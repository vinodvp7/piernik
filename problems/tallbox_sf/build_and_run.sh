#!/bin/bash
# Build and launch the tallbox_sf comparison runs (HD / MHD / diffusive CR / two-moment CR / streaming CR).
#
# Usage (from the piernik root):
#   problems/tallbox_sf/build_and_run.sh <res> <nproc> [cases...]
#     res   : 20 | 10 | 8        (cell size ~20.8 / 10.4 / 7.8 pc on the 1 x 1 x 4 kpc box)
#     nproc : MPI ranks per run
#     cases : any of A B C1 C2 D (default: all)
#   e.g.  problems/tallbox_sf/build_and_run.sh 10 64 B D
#
# Runs go to runs/tallbox_sf_<case>_<res>pc/ and are started with nohup, so they survive the terminal/session.
# Restarting after a crash or a wall-clock limit: rerun the same command; problem.par has restart = 'last'.
# Set BUILD_ONLY=1 to only compile and prepare the run directories.
#
# The dummy first define (H5PFC_EATS_1ST) works around h5pfc wrappers that drop their first argument;
# it is harmless elsewhere.

set -e

res=${1:?resolution: 20, 10 or 8}
nproc=${2:?number of MPI ranks}
shift 2
cases=${*:-A B C1 C2 D}

case $res in
   20) nd="48,48,192"   ; cred_min="1.0e4" ;;
   10) nd="96,96,384"   ; cred_min="2.0e4" ;;
   8)  nd="128,128,512" ; cred_min="2.5e4" ;;
   *)  echo "unknown resolution $res (use 20, 10 or 8)" ; exit 1 ;;
esac

declare -A defs=( [A]="H5PFC_EATS_1ST,NONMAGNETIC" [B]="H5PFC_EATS_1ST" [C1]="H5PFC_EATS_1ST,COSM_RAYS" [C2]="H5PFC_EATS_1ST,STREAM_CR" [D]="H5PFC_EATS_1ST,STREAM_CR" )

root=$(pwd)
[ -x setup ] || { echo "run this from the piernik root directory" ; exit 1 ; }

for c in $cases; do
   tag=${c}_${res}pc
   ./setup tallbox_sf -o "$tag" -d "${defs[$c]}" -p "problem.par.$c" > "build_$tag.log" 2>&1 || { echo "build $tag failed, see build_$tag.log" ; exit 1 ; }
   rundir=runs/tallbox_sf_$tag
   cp problems/tallbox_sf/cool_Kobayashi.txt "$rundir/"

   # resolution-dependent settings go through the command-line namelist override, problem.par stays untouched
   nml="&BASE_DOMAIN n_d = $nd /"
   case $c in C2|D) nml="$nml &STREAMING_CR cred_min = $cred_min /" ;; esac
   echo "$nml" > "$rundir/cmdline_nml.txt"

   if [ -z "$BUILD_ONLY" ]; then
      ( cd "$rundir" && nohup mpirun -np "$nproc" ./piernik -n "$nml" > "run_$(date +%Y%m%d_%H%M%S).log" 2>&1 & )
      echo "started $tag on $nproc ranks in $rundir"
   else
      echo "prepared $tag in $rundir (start with: cd $rundir && nohup mpirun -np $nproc ./piernik -n \"\$(cat cmdline_nml.txt)\" > run.log 2>&1 &)"
   fi
done
cd "$root"
