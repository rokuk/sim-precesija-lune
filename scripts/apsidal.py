import csv
import math

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from scipy.signal import find_peaks

import utils

# Če je True zapiše rezultate v datoteko
write = True

# Št. let
leta = 100

# Uvozi rešitve ode
intervali, rZemlja, rLuna, vZemlja, vLuna = utils.readstate(365*24*leta)

print("Računam...")

rTezisce = (rZemlja * float(utils.gmZemlja) + rLuna * float(utils.gmLuna)) \
    / (float(utils.gmZemlja) + float(utils.gmLuna))
rTezisceVektorji = np.swapaxes(rTezisce, 0, 1)
rZemljaVektorji = np.swapaxes(rZemlja, 0, 1)
rLunaVektorji = np.swapaxes(rLuna, 0, 1)
rLunaVektorjiTeziscni = [rLunaVektorji[n] - rTezisceVektorji[n] \
    for n in range(0, len(rLunaVektorji))]

# Razdalja med Zemljo in Luno
dZL = np.array([norm(rZemljaVektorji[k] - rLunaVektorji[k]) \
    for k in range(0, len(rZemljaVektorji))])[4488:]

# Poišči maksimume in minimume razdalje med Zemljo in Luno
apogeji, _ = find_peaks(dZL)
perigeji, _ = find_peaks(-1*dZL)

apogeji = np.array(apogeji) + 4488
perigeji = np.array(perigeji) + 4488

# Določi smerne vektorje apogejev in perigejev
smerniVektorjiApogeji = []
smerniVektorjiPerigeji = []

for apoUra in apogeji:
    smerniVektorjiApogeji.append(rTezisceVektorji[apoUra] - rLunaVektorji[apoUra])

for periUra in perigeji:
    smerniVektorjiPerigeji.append(rLunaVektorji[periUra] - rTezisceVektorji[periUra])

# Združi apogeje in perigeje v en array po vrstnem redu
smerniVektorji = []
maxim = len(apogeji) if len(apogeji) < len(perigeji) else len(perigeji)
for j in range(0, maxim):
    if apogeji[0] < perigeji[0]:
        smerniVektorji.append(smerniVektorjiApogeji[j])
        smerniVektorji.append(smerniVektorjiPerigeji[j])
    else:
        smerniVektorji.append(smerniVektorjiPerigeji[j])
        smerniVektorji.append(smerniVektorjiApogeji[j])

ureApsid = np.concatenate((apogeji, perigeji))
ureApsid.sort()

normala = utils.ekliptika(rTezisceVektorji[0:180*24])

# Koti med smernimi vektorji in Soncem
koti = []
tau = []
for k in range(0, len(smerniVektorji)):
    normala = np.cross(rLunaVektorjiTeziscni[ureApsid[k]-7*24], \
        rLunaVektorjiTeziscni[ureApsid[k]])
    koti.append(utils.kot(smerniVektorji[0], smerniVektorji[k], normala))
    tau.append(
        utils.kot(rTezisceVektorji[ureApsid[k]], smerniVektorji[k], normala))

ureApsid = ureApsid[0:len(koti)]

# Shrani v datoteko
if write:
    print("Writing...")
    with open(utils.path + 'data/apsidal-angles.csv', mode='w+') as outfile:
        wrt = csv.writer(outfile, delimiter=',')

        for i in range(0, len(ureApsid)):
            wrt.writerow([ureApsid[i], koti[i], tau[i]])
