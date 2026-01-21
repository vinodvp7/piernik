import sys
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 720,
        "axes.linewidth": 2.5,
        "font.size": 12,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": True,
        "legend.edgecolor": "black",
        "legend.framealpha": 1.0,
        "figure.figsize": (8, 6),
    }
)
if len(sys.argv) < 2:
    print("Usage: python plot_energy_exchange.py file.tsl")
    sys.exit(1)
df = pd.read_csv(sys.argv[1], sep=r'\s+')
t  = df['time']
ei = df['eint']
ek = df['ekin']
em = df['emag']
et = df['ener']

plt.plot(t,ei,label='eint',linewidth=1.5)
plt.plot(t,ek,label='ekin',linewidth=1.5)
plt.plot(t,em,label='emag',linewidth=1.5)
plt.plot(t,ei + ek + em,label='etot',linestyle='dashed',linewidth=1.5)
plt.plot(t,et,label='ener',linewidth=1.5)

plt.xlabel('time')
plt.ylabel('Energy')
plt.legend()
plt.savefig('energy.png',dpi=300)