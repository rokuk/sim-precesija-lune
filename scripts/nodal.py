import csv
import math

import numpy as np
from numpy.linalg import norm

import utils

# Če je True zapiše rezultate v datoteko
write = True

# Št. dni za izračun ekliptike
ekldnevi = 180

# Št. let
leta = 100

# Uvozi rešitve ode
intervali, rZemlja, rLuna, vZemlja, vLuna = utils.readstate(365*24*leta)

print("Računam...")

rTezisceZL = (float(utils.gmZemlja) * rZemlja + float(utils.gmLuna) * rLuna) \
    / (float(utils.gmZemlja) + float(utils.gmLuna))
rTezisceVektorji = np.swapaxes(rTezisceZL, 0, 1)
rTezisceVektorjiDnevi = rTezisceVektorji[0:len(rTezisceVektorji):24]
rZemljaVektorji = np.swapaxes(rZemlja, 0, 1)
rLunaVektorji = np.swapaxes(rLuna, 0, 1)
rLunaVektorjiDnevi = rLunaVektorji[0:len(rLunaVektorji):24]

# Za vsak interval po 180. dnevu od začetka do 180. dneva pred koncem
# izračunaj normalo ekliptike in skalarni produkt
# ko se predznak spremeni, dodaj na seznam dni
n = []
prehoddnevi = []
predznak = 1
for i in range(ekldnevi, 365*leta-ekldnevi+1):
    m = np.abs(np.array(utils.ekliptika(
        rTezisceVektorjiDnevi[i-ekldnevi:i+ekldnevi+1])))
    n.append(m)

    produkt = np.dot(m, rLunaVektorjiDnevi[i])

    if (produkt > 0 and predznak < 0) or (produkt < 0 and predznak > 0):
        if i != ekldnevi:
            prehoddnevi.append(i)
            predznak *= -1
        else:
            predznak *= -1

# Izračunaj prehode na uro natančno
prehodure = []
smerniVektorji = []
lunaVektorji = []
for cas in prehoddnevi:
    nh = np.linspace((n[cas-ekldnevi-2] + n[cas-ekldnevi-1])/2,
                     (n[cas-ekldnevi]+n[cas-ekldnevi+1])/2, 48)

    predznak = 1
    for ura in range(math.floor((cas-1.5)*24), math.ceil((cas+0.5)*24)):
        count = ura - math.floor((cas-1.5)*24)

        produkt = np.dot(nh[count], rLunaVektorji[ura])

        prehod = 0
        if (produkt > 0 and predznak < 0):
            prehod = 1  # gor
        elif (produkt < 0 and predznak > 0):
            if count != 0:
                prehod = 2  # dol
            else:
                predznak *= -1  # začetek

        # ko pride do prehoda
        if prehod == 1 or prehod == 2:
            prehodure.append(ura)
            lunaVektorji.append(rLunaVektorji[ura])
            predznak *= -1

        # prehod gor
        if prehod == 1:
            smerniVektorji.append(rLunaVektorji[ura] - rTezisceVektorji[ura])
        # prehod dol
        elif prehod == 2:
            smerniVektorji.append(rTezisceVektorji[ura] - rLunaVektorji[ura])

# Izračunaj kote
koti = []
kotiTezisceSmerni = []
for i in range(0, len(smerniVektorji)):
    koti.append(utils.kot(smerniVektorji[0], smerniVektorji[i], n[0]))
    kotiTezisceSmerni.append(
        utils.kot(rTezisceVektorji[prehodure[i]], smerniVektorji[i], n[0]))

# Shrani v datoteko
if write:
    print("Writing...")
    with open(utils.path + 'data/nodal-angles.csv', mode='w+') as outfile:
        wrt = csv.writer(outfile, delimiter=',')

        for i in range(0, len(koti)):
            wrt.writerow([prehodure[i], koti[i], kotiTezisceSmerni[i], \
                smerniVektorji[i][0], smerniVektorji[i][1], smerniVektorji[i][2]])
