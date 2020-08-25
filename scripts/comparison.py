import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import fftpack
import utils
from scipy.linalg import norm
from scipy.signal import find_peaks


# LaTeX fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})
plt.style.use('fast')

hZ = [[], [], []]
hL = [[], [], []]
times = []

print ("Loading Horizons data...")

with open(utils.path + 'data/horizons_results_daily100.txt','r') as fH:
    for row in fH:
        elements = row.split(',')
        times.append(float(elements[0]))
        hZ[0].append(float(elements[2]))
        hZ[1].append(float(elements[3]))
        hZ[2].append(float(elements[4]))

hZ = np.array(hZ)

intervali, srZ, srL, srvZ, srvL = utils.readstate(36525, 24)

# Preoblikuje sezname v sezname 3d vektorjev
rfZemlja = np.swapaxes(srZ, 0, 1)
rfLuna = np.swapaxes(srL, 0, 1)
vfZemlja = np.swapaxes(srvZ, 0, 1)
vfLuna = np.swapaxes(srvL, 0, 1)

# Pretvori vrednosti iz decimal v float
sZ, sL, svZ, svL = utils.makefloat(rfZemlja, rfLuna, vfZemlja, vfLuna)
sZs = np.swapaxes(sZ, 0, 1)
hZs = np.swapaxes(hZ, 0, 1)

dZh = [norm(hZs[i]-rfZemlja[i]) for i in range(0, len(hZs))]

print("Plotting...")

""" plt.figure()
plt.ylabel(u'$x_Z$ [au]')
plt.xlabel('št. dni po začetnem stanju')
plt.yticks(np.arange(-1.0, 1.1, step=0.5))
plt.plot(intervali, hZ[0], label='Horizons', marker='.', linestyle='')
plt.plot(intervali, sZs[0], label='model', marker='.', linestyle='')
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .09), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.tight_layout()
plt.savefig(utils.path + '/paper/slike/compare.eps')

plt.figure()
plt.ylabel(u'$\Delta x_Z$ [au]')
plt.xlabel('št. let po začetnem stanju')
plt.yticks([-0.005, -0.0025, 0, 0.0025, 0.005])
plt.plot(np.array(intervali)/365.2524, sZs[0]-hZ[0], color='red')
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + '/paper/slike/comparedelta.eps')
"""
plt.figure()
plt.ylabel(u'$\mid\mid \mathbf{r_Z}-\mathbf{r_{ZH}} \mid\mid$ [km]')
plt.xlabel('št. let po začetnem stanju')
plt.plot(np.array(intervali)/365.25, np.array(dZh)*float(utils.au)/1000)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + 'paper/slike/comparedistance.eps')

plt.show()
