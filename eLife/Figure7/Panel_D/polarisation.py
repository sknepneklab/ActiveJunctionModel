import numpy as np
from scipy.linalg import logm
import matplotlib.pyplot as plt
import glob
from os.path import exists

from numpy.lib import utils

from AJMAnalysis import *


def bootstrap(data, N):
    m = np.zeros(N)
    for i in range(N):
        m[i] = np.mean(np.random.choice(data, N))
    mean = np.mean(m)
    low, high = np.percentile(m, [2.5, 97.5])
    return [mean, low, high]

fpull = np.arange(0.0, 0.601, 0.025)
beta = 0.5

if not exists('fig_7D.npz'):

    Pol_xx = np.zeros((fpull.size, 3))
    Pol_yy = np.zeros((fpull.size, 3))
    Myo_xx = np.zeros((fpull.size, 3))
    Myo_yy = np.zeros((fpull.size, 3))
    T_xx = np.zeros((fpull.size, 3))
    T_yy = np.zeros((fpull.size, 3))
    for (j, fp) in enumerate(fpull):
        print(f'Processing fpull = {fp:0.4f}.')
        M0_files = sorted(
            glob.glob(f'DATA/M0_beta_{beta:0.4f}_fpull_{fp:0.4f}_seed_*.dat'))
        M_files = sorted(
            glob.glob(f'DATA/M_fpull_{fp:0.4f}_beta_{beta:0.3f}_seed_*.dat'))
        Myo_files = sorted(
            glob.glob(f'DATA/Mmyo_fpull_{fp:0.4f}_beta_{beta:0.3f}_seed_*.dat'))
        T_files = sorted(
            glob.glob(f'DATA/Tens_fpull_{fp:0.4f}_beta_{beta:0.3f}_seed_*.dat'))
        Pol_temp = np.zeros((len(M0_files), 2))
        Myo_temp = np.zeros((len(Myo_files), 2))
        T_temp = np.zeros((len(Myo_files), 2))
        for i, (m0f, mf, myof, tf) in enumerate(zip(M0_files, M_files, Myo_files, T_files)):
            M0 = np.loadtxt(m0f)[:].reshape((2,2))
            M = np.loadtxt(mf)[700, 1:].reshape((2, 2))
            Myo = np.loadtxt(myof)[700, 1:].reshape((2, 2))
            T = np.loadtxt(tf)[700, 1:].reshape((2, 2))
            U = 0.5*(logm(M) - logm(M0))
            #Pol_temp[i,:] = np.sort(np.linalg.eigh(U)[0])
            Pol_temp[i, :] = [U[0, 0], U[1, 1]]
            Myo_temp[i, :] = [Myo[0, 0], Myo[1, 1]]
            T_temp[i, :] = [T[0, 0], T[1, 1]]
        Pxx = Pol_temp[:, 0]
        Pyy = Pol_temp[:, 1]
        Pol_xx[j, :] = bootstrap(Pxx, 1000)
        Pol_yy[j, :] = bootstrap(Pyy, 1000)
        Myoxx = Myo_temp[:, 0]
        Myoyy = Myo_temp[:, 1]
        Myo_xx[j, :] = bootstrap(Myoxx, 1000)
        Myo_yy[j, :] = bootstrap(Myoyy, 1000)
        Txx = T_temp[:, 0]
        Tyy = T_temp[:, 1]
        T_xx[j, :] = bootstrap(Txx, 1000)
        T_yy[j, :] = bootstrap(Tyy, 1000)
    np.savez('fig_7D.npz', Pol_xx = Pol_xx, Pol_yy = Pol_yy, Myo_xx = Myo_xx, Myo_yy = Myo_yy, T_xx = T_xx, T_yy = T_yy) 
else: 
    data = np.load('fig_7D.npz')
    Pol_xx = data['Pol_xx']
    Pol_yy = data['Pol_yy']
    Myo_xx = data['Myo_xx']
    Myo_yy = data['Myo_yy']
    T_xx = data['T_xx']
    T_yy = data['T_yy']


plt.style.use('scientific.mplstyle')
fig, ax = plt.subplots(1, 1)
plt.errorbar(fpull, Pol_xx[:, 0], yerr=Pol_xx[:, 2] - Pol_xx[:, 1],
             linestyle='-', marker='o',  c='k', label=r'$\hat{U}_{xx}$')
plt.errorbar(fpull, Pol_yy[:, 0], yerr=Pol_yy[:, 2] - Pol_yy[:, 1],
             linestyle='--', marker='x', c='k', label=r'$\hat{U}_{yy}$')
plt.errorbar(fpull, Myo_xx[:, 0], yerr=Myo_xx[:, 2] - Myo_xx[:, 1],
             linestyle='-', marker='o', c='g', label=r'$\hat{M}^{myo}_{xx}$')
plt.errorbar(fpull, Myo_yy[:, 0], yerr=Myo_yy[:, 2] - Myo_yy[:, 1],
             linestyle='--', marker='x', c='g', label=r'$\hat{M}^{myo}_{yy}$')
plt.errorbar(fpull, T_xx[:, 0], yerr=T_xx[:, 2] - T_xx[:, 1],
             linestyle='-', marker='o', c='r', label=r'$\hat{T}_{xx}$')
plt.errorbar(fpull, T_yy[:, 0], yerr=T_yy[:, 2] - T_yy[:, 1],
             linestyle='--', marker='x', c='r', label=r'$\hat{T}_{yy}$')
ax.set_xlabel(r'$f_{pull} [f^*]$')
ax.set_ylabel('polarisation')
ax.tick_params(axis='both')
ax.set_xlim((-0.01, 0.61))
ax.set_ylim((-0.25, 1.45))
ax.set_xticks(np.arange(0.0, 0.61, 0.1))
ax.set_yticks(np.arange(-0.2, 1.41, 0.2))
plt.legend(ncol=3)
plt.savefig('fig_7D.pdf', backend='pdf')
plt.close()
