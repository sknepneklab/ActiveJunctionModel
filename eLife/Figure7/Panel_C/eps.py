import numpy as np
import matplotlib.pyplot as plt
import glob
import json
from os.path import exists
def bootstrap(data, N):
    m = np.zeros(N)
    for i in range(N):
        m[i] = np.mean(np.random.choice(data,N))
    mean = np.mean(m)
    low, high = np.percentile(m, [5,95])
    return [mean, low, high]

if not exists('vdata_to_plot.json'):
    print('Generating plot data from simulations')
    files = glob.glob('VDATA/V*.dat')

    fpull = []
    beta = []
    for f in files:
        sl = f.split('_')
        if sl[2] not in fpull:
            fpull.append(sl[2])
        if sl[4] not in beta:
            beta.append(sl[4])

    fpull.sort()
    beta.sort()

    to_plot = {}
    for b in beta:
        print(f'Processing beta = {b}.')
        to_plot[b] = {'fpull': [], 'eps': [], 'min': [], 'max': []}
        for fp in fpull:
            print(f'\tProcessing fpull = {fp}.')
            to_plot[b]['fpull'].append(float(fp))
            files = glob.glob(f'VDATA/V_fpull_{fp}_beta_{b}_*.dat')
            eps = np.zeros(len(files))
            for i,f in enumerate(files):
                data = np.loadtxt(f)
                Vxx, Vyy = np.cumsum(data[:,1]), np.cumsum(data[:,-1])
                eps[i] = Vxx[700] - Vyy[700]
            m, l, h = bootstrap(eps, 100)
            to_plot[b]['eps'].append(m)
            to_plot[b]['min'].append(l)
            to_plot[b]['max'].append(h)
else:
    with open('vdata_to_plot.json','r') as vdata:
        to_plot = json.load(vdata)
    beta = list(to_plot.keys())

plt.style.use('scientific.mplstyle')
fig, ax = plt.subplots(1,1,tight_layout=True)
for i,b in enumerate(beta): 
    plt.errorbar(to_plot[b]['fpull'], to_plot[b]['eps'], yerr = np.array(to_plot[b]['max'])-np.array(to_plot[b]['min']), marker='o', linestyle='--', label=f'{b[:-1]:s}') 
ax.set_xlabel(r'$f_{pull} [f^*]$')
ax.set_ylabel(r'$\varepsilon^{tot}_{xx}-\varepsilon^{tot}_{yy}$')
ax.set_xlim((0,0.61))
ax.set_ylim((-0.5, 0.69))
ax.set_xticks(0.1*np.arange(7))
plt.legend(ncol=4)
plt.savefig('test.pdf', backend='pdf')
plt.close()
