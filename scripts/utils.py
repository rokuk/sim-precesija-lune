import numpy as np
from numpy.linalg import norm as norm
from scipy.integrate import odeint
import csv
from decimal import *
import math

getcontext().prec = 64

"""
Initial conditions 

https://ssd.jpl.nasa.gov/horizons.cgi

JD 2458671.500000000 = A.D. 2019-Jul-07 00:00:00.0000 TDB

Zemlja
X = 2.499646121224160E-01 Y =-9.773336829448052E-01 Z = 2.463574436805279E-05
VX= 1.638044451748829E-02 VY= 4.210082689617622E-03 VZ=-4.284581571560084E-07
LT= 5.826302804458124E-03 RG= 1.008793058926798E+00 RR=-1.994876027025823E-05

Luna
X = 2.476511762731858E-01 Y =-9.765652528404914E-01 Z = 2.077343737630762E-04
VX= 1.617465022148979E-02 VY= 3.626725674604836E-03 VZ= 3.388411287849449E-05
LT= 5.818704564279764E-03 RG= 1.007477464422146E+00 RR= 4.605005402446291E-04

Coordinate system description:

  Ecliptic and Mean Equinox of Reference Epoch

    Reference epoch: J2000.0
    XY-plane: plane of the Earth's orbit at the reference epoch
              Note: obliquity of 84381.448 arcseconds wrt ICRF equator (IAU76)
    X-axis  : out along ascending node of instantaneous plane of the Earth's
              orbit and the Earth's mean equator at the reference epoch
    Z-axis  : perpendicular to the xy-plane in the directional (+ or -) sense
              of Earth's north pole at the reference epoch.

  Symbol meaning [1 au= 149597870.700 km, 1 day= 86400.0 s]:

    JDTDB    Julian Day Number, Barycentric Dynamical Time
      X      X-component of position vector (au)                               
      Y      Y-component of position vector (au)                               
      Z      Z-component of position vector (au)                               
      VX     X-component of velocity vector (au/day)                           
      VY     Y-component of velocity vector (au/day)                           
      VZ     Z-component of velocity vector (au/day)                           
      LT     One-way down-leg Newtonian light-time (day)                       
      RG     Range; distance from coordinate center (au)                       
      RR     Range-rate; radial velocity wrt coord. center (au/day)            

Geometric states/elements have no aberrations applied.
""" 

path = "/Users/rokuk/Documents/Code.nosync/Python/sim-precesija-lune/"
startdate = 2458671.5

# Parametri
au = Decimal(149597870700) # m
dan = Decimal(86400) # s

# au3 dan-2
gmSonce = Decimal(0.295912208285591100e-3)
gmZemlja = Decimal(0.888769244512563400e-9)
gmLuna = Decimal(0.109318945074237400e-10)
g = Decimal(6.6740831e-11) # m3 kg-1 s-2

rSonce = [0, 0, 0]

# Začetni pogoji, glej komentar zgoraj
# Sun (body center)
xZemlja = Decimal(2.522595272679229e-1)
yZemlja = Decimal(-9.849508544504713e-1)
zZemlja = Decimal(4.216135241692986E-5)

xLuna = Decimal(2.499460914186927e-1)
yLuna = Decimal(-9.841824243461574e-1)
zLuna = Decimal(2.252599818120031e-4)

vxZemlja = Decimal(1.638886001745495e-2)
vyZemlja = Decimal(4.210105067413899e-3)
vzZemlja = Decimal(-6.477403218368142e-7)

vxLuna = Decimal(1.618306572145645E-2)
vyLuna = Decimal(3.626748052401113E-3)
vzLuna = Decimal(3.366483071381394e-5)

# Initial state
y0 = [xZemlja, yZemlja, zZemlja, xLuna, yLuna, zLuna,\
    vxZemlja, vyZemlja, vzZemlja, vxLuna, vyLuna, vzLuna]


def vecfloat(vec):
    return np.array([float(vec[0]), float(vec[1]), float(vec[2])])


def makefloat(rrZemlja, rrLuna, vrZemlja, vrLuna):
    rfZemlja, rfLuna, vfZemlja, vfLuna = [], [], [], []

    for o in range(0, len(rrZemlja)):
        rfZemlja.append(vecfloat(rrZemlja[o]))
        rfLuna.append(vecfloat(rrLuna[o]))
        vfZemlja.append(vecfloat(vrZemlja[o]))
        vfLuna.append(vecfloat(vrLuna[o]))

    return [np.array(rfZemlja), np.array(rfLuna), np.array(vfZemlja), np.array(vfLuna)]


# Preoblikuje state v numpy vektorje, ki so bolje uporabne oblike
def repack(state):
    xZemlja, yZemlja, zZemlja, xLuna, yLuna, zLuna, \
        vxZemlja, vyZemlja, vzZemlja, vxLuna, vyLuna, vzLuna = state

    rZemlja = np.array([xZemlja, yZemlja, zZemlja])
    rLuna = np.array([xLuna, yLuna, zLuna])
    vZemlja = np.array([vxZemlja, vyZemlja, vzZemlja])
    vLuna = np.array([vxLuna, vyLuna, vzLuna])

    return [rZemlja, rLuna, vZemlja, vLuna]
    

# Vrne vsoto kinetične in potencialne energije sistema
def energija(rZemlja, rLuna, vZemlja, vLuna):
    rZemlja, rLuna, vZemlja, vLuna = \
        vecfloat(rZemlja), vecfloat(rLuna), vecfloat(vZemlja), vecfloat(vLuna)

    gmZemljaf, gmSoncef, gmLunaf, gf, auf, danf = \
        float(gmZemlja), float(gmSonce), float(gmLuna), float(g), float(au), float(dan)

    wk = 0.5*gmZemljaf/gf*auf**3/danf**2*(norm(vZemlja*auf/danf))**2 \
        + 0.5*gmLunaf/gf*auf**3/danf**2*(norm(vLuna*auf/danf))**2

    wp = - gmSoncef*gmZemljaf/gf*auf**6/danf**4/(norm(rZemlja*auf)) \
                -gmSoncef*gmLunaf/gf*auf**6/danf**4/(norm(rLuna*auf)) \
                    -gmZemljaf*gmLunaf/gf*auf**6/danf**4/(norm((rZemlja-rLuna)*auf))
    return wk + wp


# Korak
def odvajaj(t, state):
    *_r, vxZ, vyZ, vzZ, vxL, vyL, vzL = state
    rZemlja, rLuna, _vZemlja, _vLuna = repack(state)
    rZemlja = np.array([Decimal(rZemlja[0]), Decimal(rZemlja[1]), Decimal(rZemlja[2])])
    rLuna = np.array([Decimal(rLuna[0]), Decimal(rLuna[1]), Decimal(rLuna[2])])
    
    zl = rZemlja - rLuna
    nd = Decimal(norm(zl))
    nz = Decimal(norm(rZemlja))
    nl = Decimal(norm(rLuna))
    dolzinaZL = nd*nd*nd

    # Vektorski račun!
    aZ = -gmSonce*rZemlja/(nz*nz*nz) - gmLuna*zl/dolzinaZL
    aL = -gmSonce*rLuna/(nl*nl*nl) + gmZemlja*zl/dolzinaZL

    return [vxZ, vyZ, vzZ, vxL, vyL, vzL, aZ[0], aZ[1], aZ[2], aL[0], aL[1], aL[2]]


# Reši ode
def solve(start, stop, N):
    print("Solving...")

    # Časi ob katerih izračunamo rešitve
    intervali = np.linspace(start, stop, num=N)

    # Reši enačbe
    result = odeint(odvajaj, y0, intervali, tfirst=True, ixpr=True, \
        full_output=False, hmax=0.02, hmin=0.01, rtol=1.0e-32)

    # Seznami rešitev
    return [intervali, *repack([result[:,0], result[:,1], \
        result[:,2], result[:,3], result[:,4], result[:,5], result[:,6], \
            result[:,7], result[:,8], result[:,9], result[:,10], result[:,11]])]


# Izračuna vektorje pri podanih časovnih intervalih - POZOR! Manjša natančnost!
def solveinterval(intervali):
    print("Solving...")

    # Reši enačbe
    result = odeint(odvajaj, y0, intervali, tfirst=True, ixpr=True, \
        full_output=False, hmax=0.05, rtol=1.0e-32)

    # Seznami rešitev
    return [intervali, *repack([result[:,0], result[:,1], \
        result[:,2], result[:,3], result[:,4], result[:,5], result[:,6], \
            result[:,7], result[:,8], result[:,9], result[:,10], result[:,11]])]


# Izračuna seznam normal ekliptike
def ekliptika(rZemljaVektorji):
    G = rZemljaVektorji.sum(axis=0) / rZemljaVektorji.shape[0] # težišče točk
    u, s, vh = np.linalg.svd(rZemljaVektorji - G)
    return vh[2, :] # normala ekliptike


# Izračuna kot med vektorjema v kodomeni 0 in 2pi okrog normale
def kot(e1, en, n):
    number = np.dot(e1, en)/(norm(e1)*norm(en))
    if number > 1: number = 1
    acosdot = math.acos(number)
    smer = np.dot(np.cross(e1, en), n)
    kot = 0

    if smer > 0:
        kot = acosdot
    elif smer < 0:
        kot = 2*math.pi - acosdot
    
    return kot


# Shrani state vektorje v datoteko
def savestate(state):
    print("Writing state...")

    intervali, rZemlja, vZemlja, rLuna, vLuna = state

    with open(path + 'data/state.csv', mode='w+') as outfile:
        wrt = csv.writer(outfile, delimiter=',')

        for i in range(0, len(intervali)):
            wrt.writerow([intervali[i], \
                rZemlja[0][i], rZemlja[1][i], rZemlja[2][i],\
                    vZemlja[0][i], vZemlja[1][i], vZemlja[2][i],\
                        rLuna[0][i], rLuna[1][i], rLuna[2][i],\
                            vLuna[0][i], vLuna[1][i], vLuna[2][i]])


# Uvozi state vektorje iz datoteke
def readstate(n, step=1):
    print("Importing state...")

    intervali = []
    rZemlja = [[], [], []]
    vZemlja = [[], [], []]
    rLuna = [[], [], []]
    vLuna = [[], [], []]
    
    with open(path + 'data/state.csv') as csv_file:
        reader = csv.reader(csv_file, delimiter=',')

        for i in range(0, n*step):
            row = reader.__next__()
            if i % step == 0:
                intervali.append(float(row[0]))
                rZemlja[0].append(float(row[1]))
                rZemlja[1].append(float(row[2]))
                rZemlja[2].append(float(row[3]))
                vZemlja[0].append(float(row[4]))
                vZemlja[1].append(float(row[5]))
                vZemlja[2].append(float(row[6]))
                rLuna[0].append(float(row[7]))
                rLuna[1].append(float(row[8]))
                rLuna[2].append(float(row[9]))
                vLuna[0].append(float(row[10]))
                vLuna[1].append(float(row[11]))
                vLuna[2].append(float(row[12]))

    return [np.array(intervali), np.array(rZemlja), np.array(vZemlja), \
        np.array(rLuna), np.array(vLuna)]
