#!/usr/bin/env python3
import sys, os, argparse
import h5py, numpy as np
import xml.etree.ElementTree as ET

def kji_str(dims):
    if len(dims)==3:
        nx,ny,nz = dims
        return f"{nz} {ny} {nx}"
    if len(dims)==2:
        nx,ny = dims
        return f"{ny} {nx}"
    if len(dims)==1:
        return f"{dims[0]}"
    raise ValueError("bad dims")

def topo_type(nd): return "3DCORECTMesh" if nd==3 else ("2DCORECTMesh" if nd==2 else None)
def geom_type(nd): return "ORIGIN_DXDYDZ" if nd==3 else ("ORIGIN_DXDY" if nd==2 else None)

def numtype(dt):
    if dt.kind=="f": return "Float"
    if dt.kind in ("i","u"): return "Int"
    return "Float"
def prec(dt):
    if dt.kind=="f" and dt.itemsize==8: return "8"
    if dt.kind=="f" and dt.itemsize==4: return "4"
    return None

def prod(a):
    p=1
    for v in a: p*=int(v)
    return p

def infer_ng(shape, nb):
    """Symmetric ghost detection; returns ng per axis or None."""
    if len(shape)!=len(nb): return None
    ng=[]
    for s,c in zip(shape, nb):
        d=s-c
        if d<0 or d%2!=0: return None
        ng.append(d//2)
    return ng

def choose_squeeze_axis(nb, mode):
    if mode=="none" or len(nb)==2: return None
    if mode in ("x","y","z"): return {"x":0,"y":1,"z":2}[mode]
    cand=[i for i,n in enumerate(nb) if n==1]
    return cand[-1] if cand else None

def build_collection(root, name, gridspecs, topology, log):
    coll = ET.SubElement(root,"Grid", Name=name, GridType="Collection", CollectionType="Spatial")
    for spec in gridspecs:
        # gridspec can be (gname, nb, le, dl, vars) or (gname, nb, le, dl, vars, meta)
        if len(spec) == 5:
            (gname, nb, le, dl, vars_out) = spec
        else:
            (gname, nb, le, dl, vars_out, meta) = spec
        nd=len(nb)
        topo=topo_type(nd); geom=geom_type(nd)
        gxml=ET.SubElement(coll,"Grid", Name=gname, GridType="Uniform")
        if topology=="nodes":
            nodes=[n+1 for n in nb]
            ET.SubElement(gxml,"Topology", Name=name, TopologyType=topo, Dimensions=kji_str(nodes))
        else:
            ET.SubElement(gxml,"Topology", Name=name, TopologyType=topo, Dimensions=kji_str(nb))
        gmx=ET.SubElement(gxml,"Geometry", Name=name, GeometryType=geom)
        ET.SubElement(gmx,"DataItem", Dimensions=str(nd), NumberType="Float", Format="XML").text=" ".join(map(str, le[:nd]))
        ET.SubElement(gmx,"DataItem", Dimensions=str(nd), NumberType="Float", Format="XML").text=" ".join(map(str, dl[:nd]))

        for v in vars_out:
            dname, src_dims_kji, target_dims_kji, src_path, nt, pr, slab_meta = v
            attr=ET.SubElement(gxml,"Attribute", Name=dname, AttributeType="Scalar", Center="Cell")
            if slab_meta is None:
                di_kwargs={"Dimensions":target_dims_kji, "NumberType":nt, "Format":"HDF"}
                if pr: di_kwargs["Precision"]=pr
                di=ET.SubElement(attr,"DataItem", **di_kwargs)
                di.text=src_path
            else:
                (slab_dims, slab_vals)=slab_meta
                di_hs=ET.SubElement(attr,"DataItem", ItemType="HyperSlab", Dimensions=target_dims_kji, Type="HyperSlab")
                ET.SubElement(di_hs,"DataItem", Dimensions=slab_dims, NumberType="Int", Format="XML").text = " ".join(map(str, slab_vals))
                src_kwargs={"Dimensions":src_dims_kji, "Format":"HDF", "NumberType":nt}
                if pr: src_kwargs["Precision"]=pr
                src=ET.SubElement(di_hs,"DataItem", **src_kwargs)
                src.text=src_path
    return coll


def rank_from_attrs(attrs):
    keys = ("mpi_rank","rank","proc","owner","pe","pe_owner")
    for k in keys:
        if k in attrs:
            try:
                v = attrs[k]
                return int(v[0] if hasattr(v,"__len__") else v)
            except Exception:
                pass
    return -1

def build_block_outlines(root, gridspecs, topology, log):
    # Create a collection of unstructured line meshes (one per block)
    coll = ET.SubElement(root, "Grid", Name="BlockOutlines", GridType="Collection", CollectionType="Spatial")
    for spec in gridspecs:
        (gname, nb, le, dl, vars, meta) = spec
        nd = len(nb)
        if nd != 2:
            continue  # outlines are 2D only
        # Compute four corners
        nx, ny = nb
        x0, y0 = le[0], le[1]
        x1, y1 = x0 + nx*dl[0], y0 + ny*dl[1]
        # Points array (4 x 2): (x0,y0),(x1,y0),(x1,y1),(x0,y1)
        pts = [x0, y0,  x1, y0,  x1, y1,  x0, y1]
        # Connectivity for 4 edges (Edge_2): [0,1],[1,2],[2,3],[3,0]
        conn = [0,1, 1,2, 2,3, 3,0]
        # Rank attribute if known
        rank = meta.get("rank", -1)

        gxml = ET.SubElement(coll, "Grid", Name=f"{gname}_outline", GridType="Uniform")
        topo = ET.SubElement(gxml, "Topology", TopologyType="Edge_2", NumberOfElements="4")

        # Connectivity as DataItem (Nelements x 2) in XML (small)
        ET.SubElement(topo, "DataItem", Dimensions="4 2", NumberType="Int", Format="XML").text = " ".join(map(str, conn))

        geom = ET.SubElement(gxml, "Geometry", GeometryType="XYZ")
        pts3 = []
        for i in range(0, len(pts), 2):
            pts3.extend([pts[i], pts[i+1], 0.0])
        ET.SubElement(geom, "DataItem", Dimensions="4 3", NumberType="Float", Format="XML").text = " ".join(map(str, pts3))

        # Optional per-edge rank (repeat same rank for 4 edges)
        attr = ET.SubElement(gxml, "Attribute", Name="rank", AttributeType="Scalar", Center="Cell")
        ET.SubElement(attr, "DataItem", Dimensions="4", NumberType="Int", Format="XML").text = " ".join([str(rank)]*4)

    return coll


def make(h5path, outpath=None, topology="nodes", squeeze="auto", emit_blocks=True, log=True):
    if outpath is None: outpath = os.path.splitext(h5path)[0] + ".xdmf"
    logpath = os.path.splitext(outpath)[0] + ".log" if log else None

    with h5py.File(h5path,"r") as f, open(logpath,"w") if log else open(os.devnull,"w") as L:
        if "data" not in f: raise RuntimeError("No '/data' group")
        gdata=f["data"]
        grids=[k for k,v in gdata.items() if isinstance(v, h5py.Group) and k.startswith("grid_")]
        grids.sort()

        time_val=None
        if "simulation_parameters" in f:
            gsp=f["simulation_parameters"]
            if "current_time" in gsp.attrs:
                tv=gsp.attrs["current_time"]
                try: time_val=float(tv[0] if hasattr(tv,"__len__") else tv)
                except: pass

        # Gather per-grid specs once; we will emit twice (AMR and Blocks) if requested
        gridspecs=[]
        for gname in grids:
            gg=gdata[gname]
            if not all(a in gg.attrs for a in ("n_b","left_edge","dl")): 
                print(f"[skip] {gname}: missing one of n_b,left_edge,dl", file=L); 
                continue

            nb = list(map(int,   np.array(gg.attrs["n_b"]).tolist()))
            le = list(map(float, np.array(gg.attrs["left_edge"]).tolist()))
            dl = list(map(float, np.array(gg.attrs["dl"]).tolist()))
            axis_squeeze = choose_squeeze_axis(nb, squeeze)
            if axis_squeeze is not None and len(nb)==3:
                nb2=[nb[i] for i in range(3) if i!=axis_squeeze]
                le2=[le[i] for i in range(3) if i!=axis_squeeze]
                dl2=[dl[i] for i in range(3) if i!=axis_squeeze]
            else:
                nb2, le2, dl2 = nb, le, dl

            cell_prod = prod(nb)
            vars_out=[]
            meta={'rank': rank_from_attrs(gg.attrs)}
            for dname, obj in gg.items():
                if not isinstance(obj, h5py.Dataset): continue
                shp=list(obj.shape); shp_kji=shp[::-1]
                nt=numtype(obj.dtype); pr=prec(obj.dtype)
                direct=False; crop=False

                if prod(shp)==cell_prod:
                    # direct mapping
                    direct=True
                else:
                    ng = infer_ng(shp, nb) or infer_ng(shp[::-1], nb)
                    if ng is not None: crop=True

                target_dims_kji = kji_str(nb2)

                if direct:
                    src_path=f"{os.path.basename(h5path)}:/data/{gname}/{dname}"
                    vars_out.append((dname, kji_str(shp_kji), target_dims_kji, src_path, nt, pr, None))
                else:
                    if infer_ng(shp, nb): ng_use = infer_ng(shp, nb)
                    else:
                        ng_use = infer_ng(shp[::-1], nb)
                        shp = shp[::-1]; shp_kji = shp_kji[::-1]
                    src_nd = len(shp_kji)
                    if src_nd==3:
                        start=[ng_use[2], ng_use[1], ng_use[0]]
                        count=[nb[2], nb[1], nb[0]]
                        slab_dims="3 3"; slab_vals=start+[1,1,1]+count
                    else:
                        start=[ng_use[1], ng_use[0]]
                        count=[nb[1], nb[0]]
                        slab_dims="3 2"; slab_vals=start+[1,1]+count
                    src_path=f"{os.path.basename(h5path)}:/data/{gname}/{dname}"
                    vars_out.append((dname, kji_str(shp_kji), target_dims_kji, src_path, nt, pr, (slab_dims, slab_vals)))
            gridspecs.append((gname, nb2, le2, dl2, vars_out, meta))

        # Build XML
        xdmf=ET.Element("Xdmf", Version="3.0")
        dom=ET.SubElement(xdmf,"Domain")
        if time_val is not None:
            # Add time info to each collection via a Time child
            coll_amr = build_collection(dom, "AMR", gridspecs, topology, L)
            ET.SubElement(coll_amr,"Time", Value=str(time_val))
            if emit_blocks:
                coll_blk = build_collection(dom, "Blocks", gridspecs, topology, L)
                ET.SubElement(coll_blk,"Time", Value=str(time_val))
            # Always include outlines for easier inspection
            coll_ol = build_block_outlines(dom, gridspecs, topology, L)
            ET.SubElement(coll_ol, "Time", Value=str(time_val))
        else:
            build_collection(dom, "AMR", gridspecs, topology, L)
            if emit_blocks: build_collection(dom, "Blocks", gridspecs, topology, L)
            build_block_outlines(dom, gridspecs, topology, L)

        # pretty print
        def indent(e, lvl=0):
            i="\n"+lvl*"  "
            if len(e):
                if not e.text or not e.text.strip(): e.text=i+"  "
                for c in e: indent(c, lvl+1)
                if not c.tail or not c.tail.strip(): c.tail=i
            if lvl and (not e.tail or not e.tail.strip()): e.tail=i
        indent(xdmf)
        ET.ElementTree(xdmf).write(outpath, encoding="utf-8", xml_declaration=True)
        return outpath, logpath

if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("xdmf", nargs="?")
    ap.add_argument("--topology", choices=("nodes","cells"), default="nodes")
    ap.add_argument("--squeeze", choices=("auto","none","x","y","z"), default="auto")
    ap.add_argument("--no-blocks", action="store_true", help="Do not emit the duplicate 'Blocks' mesh")
    args=ap.parse_args()
    out, log = make(args.h5, args.xdmf, topology=args.topology, squeeze=args.squeeze, emit_blocks=(not args.no_blocks))
    print("Wrote:", out)
    print("Log:  ", log)
