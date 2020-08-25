import csv
import math

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from scipy.signal import detrend, find_peaks

import utils


# LaTeX fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


print("Importing...")

koti = []
ure = []
tau = []
with open(utils.path + 'data/apsidal-angles.csv') as csv_file:
    reader = csv.reader(csv_file, delimiter=',')

    for row in reader:
        ure.append(int(row[0]))
        koti.append(float(row[1]))
        tau.append(float(row[2]))

leta = np.array(ure)/24/365.25

print("Plotting...")
""" 
plt.figure(figsize=(11,4.9))
plt.ylabel(r'$\gamma_i$')
plt.xlabel('št. let po začetnem stanju')
plt.yticks([0, math.pi/2, math.pi, math.pi/2*3, math.pi*2], \
    [0, r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
plt.plot(leta, koti)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + 'paper/slike/apsidal-angle.eps')

# Autocorrelation
plt.figure()
plt.xlim((0, 2600))
plt.xticks([0, 500, 1000, 1500, 2000, 2500])
plt.xlabel(r'$\Delta i$')
plt.ylabel(r'avtokorelacija $\gamma_i$')
lags, c, line, b = plt.acorr(koti, usevlines=False, marker='.', maxlags=None)
line.set_markersize(1)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/apsidal-angle-acorr.eps")

ure = np.array(ure)
#print([0, 236, 471, 704, 939, 1174, 1409, 1644, 1878, 2112, 2347])
#print(ure[[0, 236, 471, 704, 939, 1174, 1409, 1644, 1878, 2112, 2347]])

# Detrend
plt.figure(figsize=(10.6, 4))
plt.ylabel(r'$\gamma_i - \textrm{trend}$')
plt.xlabel('št. let po začetnem stanju')
plt.yticks([-0.5, -0.25, 0, 0.25, 0.5])
plt.plot(leta[0:200:2], detrend(koti[0:200:2]), marker='.', linewidth=0.9, linestyle=':')
plt.plot(leta[1:200:2], detrend(koti[1:200:2]), marker='.', linewidth=0.9, linestyle=':', color='tab:red')
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/apsidal-angle-detrend.eps")

# Detrend autocorrelation
plt.figure(figsize=(10.6, 4.2))
plt.xlim((0, 200))
plt.ylim((-1.1, 1.1))
plt.yticks([-1, -0.5, 0, 0.5, 1])
plt.xlabel(r'$\Delta i$')
plt.ylabel(r'avtokorelacija ($\gamma_i - \mathrm{trend}$)', labelpad=5)
lags, c, line, b = plt.acorr(detrend(koti[0:200]), maxlags=None)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/apsidal-angle-detrend-acorr.eps")
 """
# Primerjava
fig, (ax2, ax4) = plt.subplots(2, 1, sharex=True, figsize=(9.3,6.4))

ax2.set_ylabel(r'$\tau_i$', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.plot(leta[244:336], tau[244:336], color='tab:red', linestyle=':', marker='.')
ax2.set_yticks([0, math.pi/2, math.pi, math.pi/2*3, math.pi*2])
ax2.set_yticklabels([0, r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax2.grid(axis='y')

ax1 = ax2.twinx()
ax1.set_ylabel(r'$\gamma_i - \textrm{trend}$', color='#00578b')
ax1.tick_params(axis='y', labelcolor='#00578b')
ax1.set_ylim((-0.6, 0.6))
ax1.plot(leta[244:336], detrend(koti[244:336]), linestyle=':', marker='.')

ax4.set_ylabel(r'$\tau_i$', color='tab:red')
ax4.tick_params(axis='y', labelcolor='tab:red')
ax4.plot(leta[244:336], tau[244:336], color='tab:red', linestyle=':', marker='.')
ax4.set_yticks([0, math.pi/2, math.pi, math.pi/2*3, math.pi*2])
ax4.set_yticklabels([0, r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax4.grid(axis='y')
ax4.set_xlabel('št. let po začetnem stanju')

ax3 = ax4.twinx()
ax3.set_ylabel(r'$\gamma_i (\mathrm{perigej})$', color='#00578b')
ax3.tick_params(axis='y', labelcolor='#00578b')
ax3.set_ylim((0.4,3.5))
ax3.plot(leta[244:336:2], koti[244:336:2], linestyle=':', marker='.')

ax5 = ax4.twinx()
ax5.set_ylim((0.4,3.5))
ax5.set_ylabel(r'$\gamma_i (\mathrm{apogej})$', color='tab:orange')
ax5.tick_params(axis='y', labelcolor='tab:orange')
ax5.plot(leta[245:336:2], koti[245:336:2], linestyle=':', marker='.', color='tab:orange')

xcoords = [9.74, 10.03, 10.31, 10.59, 10.87, 11.16, 11.46, 11.73, 12.01, 12.28, 12.57, 12.84, 13.14]
for xc in xcoords:
    ax4.axvline(x=xc, color='#c9c9c9', linewidth=0.3)
    ax2.axvline(x=xc, color='#c9c9c9', linewidth=0.3)

ax5.spines["right"].set_position(("axes", 1.13))
make_patch_spines_invisible(ax5)
ax5.spines["right"].set_visible(True)

fig.savefig(utils.path + "paper/slike/apsidal-angle-comp.eps")

plt.show()
