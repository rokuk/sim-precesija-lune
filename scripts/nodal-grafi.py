import math
import csv

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from scipy.signal import detrend, find_peaks

import utils


# LaTeX fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})


print("Importing...")

koti = []
prehodure = []
kotiTezisceSmerni = []
with open(utils.path + 'data/nodal-angles.csv') as csv_file:
    reader = csv.reader(csv_file, delimiter=',')

    for row in reader:
        koti.append(float(row[1]))
        prehodure.append(int(row[0]))
        kotiTezisceSmerni.append(float(row[2]))

print("Plotting...")

plt.figure(figsize=(11,4.5))
plt.ylabel(r'$\phi_i$')
plt.xlabel('št. let po začetnem stanju')
plt.yticks([0, math.pi/2, math.pi, math.pi/2*3, math.pi*2], \
    [0, r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
plt.plot(np.array(prehodure)/24/365.25, koti)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/nodal-angle.eps")

# Autocorrelation
plt.figure()
plt.xlim((0, 2600))
plt.xticks([0, 500, 1000, 1500, 2000, 2500])
plt.xlabel(r'$\Delta i$')
plt.ylabel(r'avtokorelacija $\phi_i$')
lags, c, line, b = plt.acorr(koti, usevlines=False, marker='.', maxlags=None)
line.set_markersize(1)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/nodal-angle-acorr.eps")

# Detrend
plt.figure(figsize=(10.6, 4))
plt.ylabel(r'$\phi_i - \textrm{trend}$')
plt.xlabel('št. let po začetnem stanju')
plt.yticks([-0.02, 0, 0.02])
plt.plot(np.array(prehodure[10:-10])/24/365.25, detrend(koti[10:-10]), marker='.', linewidth=0.9)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/nodal-angle-detrend.eps")

plt.figure(figsize=(10.3, 4.5))
plt.xlabel(r'$\Delta i$')
plt.ylabel(r'avtokorelacija ($\phi_i - \textrm{trend}$)')
plt.yticks([-1, -0.5, 0, 0.5, 1])
lags, c, line, b = plt.acorr(detrend(koti[10:-10]), maxlags=None)
plt.xlim(0, len(lags)//2)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/nodal-angle-detrend-acorr.eps")

fig, (ax2, ax4) = plt.subplots(2, 1, sharex=True, figsize=(9.3,6.4))

ax2.set_ylabel(r'$\beta_i$', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.plot(np.array(prehodure[100:-330])/24/365.25, kotiTezisceSmerni[100:-330], color='tab:red', linestyle=':', marker='.')
ax2.set_yticks([0, math.pi/2, math.pi, math.pi/2*3, math.pi*2])
ax2.set_yticklabels([0, r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax2.grid(axis='y')

ax1 = ax2.twinx()
ax1.set_ylabel(r'$\phi_i - \textrm{trend}$', color='#00578b', labelpad=10)
ax1.tick_params(axis='y', labelcolor='#00578b')
ax1.plot(np.array(prehodure[100:-330])/24/365.25, detrend(koti[100:-330]), linestyle=':', marker='.')

ax4.set_ylabel(r'$\beta_i$', color='tab:red')
ax4.tick_params(axis='y', labelcolor='tab:red')
ax4.plot(np.array(prehodure[100:-330])/24/365.25, kotiTezisceSmerni[100:-330], color='tab:red', linestyle=':', marker='.')
ax4.set_yticks([0, math.pi/2, math.pi, math.pi/2*3, math.pi*2])
ax4.set_yticklabels([0, r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax4.grid(axis='y')
ax4.set_xlabel('št. let po začetnem stanju')

ax3 = ax4.twinx()
ax3.set_ylabel(r'$\phi_i$', color='#00578b', labelpad=10)
ax3.tick_params(axis='y', labelcolor='#00578b')
ax3.plot(np.array(prehodure[100:-330])/24/365.25, koti[100:-330], linestyle='', marker='.')

xcoords = [4.274, 4.497, 4.757, 4.98, 5.24, 5.47, 5.7, 5.93, 6.16, 6.432, 6.64, 6.89, 7.09]
for xc in xcoords:
    ax4.axvline(x=xc, color='#c9c9c9', linewidth=0.3, alpha=0.3)
    ax2.axvline(x=xc, color='#c9c9c9', linewidth=0.3, alpha=0.3)

fig.tight_layout()
fig.savefig(utils.path + "paper/slikep/nodal-angle-comp.eps")

#plt.show()
