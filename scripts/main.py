import math

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from scipy.signal import find_peaks

import utils

# LaTeX fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

# Prikaz grafov in shranjevanje izračunanih vektorjev v datoteko: 1: da, 0: ne
plotting = 1
fileoutput = 0

# Reši ode
intervali, rZemlja, rLuna, vZemlja, vLuna = utils.readstate(365*24*101)

# Preoblikuj sezname v sezname 3d vektorjev
rZemljaVektorji = np.swapaxes(rZemlja, 0, 1)
vZemljaVektorji = np.swapaxes(vZemlja, 0, 1)
rLunaVektorji = np.swapaxes(rLuna, 0, 1)
vLunaVektorji = np.swapaxes(vLuna, 0, 1)

# Izračunaj razdaljo med Zemljo in Luno
dZLdelt = rZemlja - rLuna
dZemljaLuna = [norm([dZLdelt[0][i], dZLdelt[1][i], dZLdelt[2][i]]) \
    for i in range(0, len(dZLdelt[0]))]

# Izračunaj težišče sistema Zemlja-Luna
rTezisceZL = (float(utils.gmZemlja) * rZemlja + float(utils.gmLuna) * rLuna) \
    / (float(utils.gmZemlja) + float(utils.gmLuna))

# Izračunaj krajevni vektor Lune in Zemlje glede na koordinatni sistem s centrom v težišču sistema Zemlja-Luna
rLunaTeziscni = rLuna - rTezisceZL
rZemljaTeziscni = rZemlja - rTezisceZL

# Obhodni časi
peaksLuna, _ = find_peaks(rLunaTeziscni[0])
peaksZemlja, _ = find_peaks(rZemlja[0])
print((np.array(intervali))[peaksZemlja])
print(len(peaksZemlja))

""" # Izračunaj energije za vsak korak
energ = np.array([utils.energija(\
    rZemljaVektorji[j], rLunaVektorji[j], vZemljaVektorji[j], vLunaVektorji[j]) \
        for j in range(0, len(rZemljaVektorji))])
deltaenerg = energ - energ[0]

# Velikost krajevnega vektorja Zemlje
drZemlja = [norm(rZemljaVektorji[k]) for k in range(0, len(rZemljaVektorji))]
 """
print("Plotting...")

"""
# Obhodni časi Lune
peakpointsLuna = []
numsLuna = []
for i in range(1, len(peaksLuna)):
    numsLuna.append(i)
    peakpointsLuna.append(intervali[peaksLuna[i]]- intervali[peaksLuna[i-1]])

plt.figure()
plt.xlabel('število periode')
plt.ylabel('dolžina Lunine periode [dni]')
plt.scatter(numsLuna, peakpointsLuna, s=2)
plt.tight_layout()
plt.savefig(utils.path + 'paper/slike/obhodi-luna.eps')

# Obhodni časi Zemlje
peakpointsZemlja = []
numsZemlja = []
for i in range(1, len(peaksZemlja)):
    numsZemlja.append(i)
    peakpointsZemlja.append(intervali[peaksZemlja[i]]- intervali[peaksZemlja[i-1]])

plt.figure()
plt.xlabel('število periode')
plt.ylabel('dolžina Zemljine periode [dni]')
plt.ylim((365.2, 365.35))
plt.yticks([365.2, 365.25, 365.3, 365.35])
plt.scatter(numsZemlja, peakpointsZemlja, s=2)
plt.tight_layout()
plt.savefig(utils.path + 'paper/slike/obhodi-zemlja.eps')

# Grafi
fig, ax = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 5))
ax[0, 0].set_ylabel(r'$x_Z$ [au]')
ax[0, 0].plot(intervali, rZemlja[0])
ax[0, 1].yaxis.tick_right()
ax[0, 1].yaxis.set_ticks_position('both')
ax[0, 1].yaxis.set_label_coords(1.25,0.5)
ax[0, 1].set_ylabel(r'$y_Z$ [au]')
ax[0, 1].plot(intervali, rZemlja[1])
ax[1, 0].set_ylabel(r'$x_{LT}$ [km]')
ax[1, 0].plot(intervali, rLunaTeziscni[0]*float(utils.au)/1000)
ax[1, 1].yaxis.tick_right()
ax[1, 1].yaxis.set_ticks_position('both')
ax[1, 1].set_ylabel(r'$y_{LT}$ [km]')
ax[1, 1].plot(intervali, rLunaTeziscni[1]*float(utils.au)/1000)
ax[1, 1].yaxis.set_label_coords(1.25,0.5)
fig.text(0.51, 0.04, 'št. dni po začetnem stanju', ha='center', va='center')
ax[1, 0].xaxis.set_ticks([0, 100, 200, 300])
ax[1, 1].xaxis.set_ticks([0, 100, 200, 300])
plt.tight_layout()
fig.subplots_adjust(hspace=0.25, wspace=0.15)
plt.savefig(utils.path + "paper/slike/koordinate.eps")

# Koordinata z
fig, ax1 = plt.subplots(figsize=(10, 3.5))

ax1.set_ylabel(r'$z_Z$ [km]')
ax1.plot(intervali, rZemlja[2]*float(utils.au)/1000, color='black')

ax2 = ax1.twinx()
ax2.set_ylabel(r'$z_L$ [km]', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.plot(intervali, rLuna[2]*float(utils.au)/1000, color='tab:red')
ax1.set_xlabel('št. dni po začetnem stanju')

fig.tight_layout()
fig.savefig(utils.path + 'paper/slike/koordinatez.eps')

plt.figure()
plt.plot(intervali/365.25, np.array(dZemljaLuna)*float(utils.au)/1000)
plt.ylabel(r'$\mid \mathbf{r_Z}-\mathbf{r_L} \mid$ [km]')
plt.xlabel('št. let po začetnem stanju')
plt.xticks([0, 25, 50 , 75, 100])
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/razdaljaZL.eps")

plt.figure()
plt.plot(intervali/365.25, np.array(drZemlja))
plt.ylabel(r'$\mid \mathbf{r_Z} \mid$ [au]')
plt.xlabel('št. let po začetnem stanju')
plt.xticks([0, 25, 50 , 75, 100])
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/razdaljaSZ.eps")

plt.figure()
plt.plot(rZemljaTeziscni[0], rZemljaTeziscni[1])
plt.plot(rLunaTeziscni[0], rLunaTeziscni[1])
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/koord-teziscni.eps")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(rZemlja[0], rZemlja[1], rZemlja[2])
ax.plot(rLuna[0], rLuna[1], rLuna[2])
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/3d-koord.eps")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(rZemljaTeziscni[0], rZemljaTeziscni[1], rZemljaTeziscni[2])
ax.plot(rLunaTeziscni[0], rLunaTeziscni[1], rLunaTeziscni[2])
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/3d-koord-teziscni.eps")
#ax.set(xlim=(-0.003, 0.003), ylim=(-0.003, 0.003), zlim=(-0.003, 0.003))

plt.figure()
plt.ylabel(u'$W - W_0$ [J]')
plt.xlabel('št. let po začetnem stanju')
plt.plot(intervali/365.25, deltaenerg)
plt.grid()
plt.tight_layout()
plt.savefig(utils.path + "paper/slike/energija.eps")
"""
print("Done!")

#print("začetna energija:", energ[0])

if plotting: plt.show()
