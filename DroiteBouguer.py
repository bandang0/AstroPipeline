#!/usr/bin/env python3

import math
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds
import PythonPhot as pp
import numpy as np

from general import *

#gaussB, psfB, psfmagB = pp.getpsf (imageB, xstarSumBV [ipal_B],  ypos,nmag,sky,1,1,np.arange(len(xpos)),5,'output_psf.fits')

def flux(name):
    print(name)
    image = fits.getdata (name)
    
    # On fit l'histogramme du fond de ciel avec une gaussienne
    skymod, skysig, skyskw = pp.mmm.mmm (image, readnoise = 0, minsky = 200)
    
    # On regarde à 5 sigma au dessus du fond de ciel pour identifier les étoiles
    hmin = skysig * 6
    
    # Largeur à mi-hauteur totale d'une étoile
    fwhm = 12
    
    # On repère les étoiles dans les images
    xstar, ystar, cflux, sharp, round = pp.find.find (image, hmin, fwhm, [-1.0,1.0], [0.2,1.0])
   
    gain = 1.0
    
    # On utilise la fonction aper pour récupérer les magnitudes et les valeurs du fond de ciel pour les étoiles repérées
    mag, magerr, cflux, fluxerr, sky, skyerr, badflag, outstr = \
            pp.aper.aper (image / gain, xstar , ystar, phpadu = gain, apr = 20, zeropoint = 25.11, skyrad = [3 * fwhm, 5 * fwhm], badpix = [-12000, 60000], exact = True, setskyval = 0)
    print(cflux)
    return float(max(cflux))

fluxesB = list()
fluxesV = list()

starsB = [photomDir + i for i in ["HD147394-B-1.fits",
                                  "HD147394-B-2.fits",
                                  "HD147394-B-3.fits"]]
starsV = [photomDir + i for i in ["HD147394-V-1.fits",
                                  "HD147394-V-2.fits",
                                  "HD147394-V-3.fits"]]
for name in starsB:
    l = flux(name)
    fluxesB.append(l)
    
for name in starsV:
    l = flux(name)
    fluxesV.append(l)
    
fluxesB = [float(i) for i in fluxesB]

fluxesV = [float(i) for i in fluxesV]


dt = [0.5, 0.5, 1.0]

magsB = [mag(f, j) for (f, j) in zip(fluxesB, dt)]
#magsBerrs = [math.fabs(mag(f + df, j) - mag(f, j)) for (f, df, j) in zip(fluxesB, fluxBerrs, dt)]

magsV = [mag(f, j) for (f, j) in zip(fluxesV, dt)]
#magsVerrs = [math.fabs(mag(f + df, j) - mag(f, j)) for (f, df, j) in zip(fluxesV, fluxVerrs, dt)]

X = [1.0/math.cos(rad(d,m)) for (d, m) in zip([57, 62, 67], [56, 10, 30])]# 67], [56, 10, 30])]


plt.plot(X, magsB, 'r+')
plt.xlabel('$1/\cos(z)$')
plt.ylabel('$Mag_B$')
plt.show()


plt.plot(X, magsV, 'r+')
plt.xlabel('$1/\cos(z)$')
plt.ylabel('$Mag_V$')
plt.text("hello")
plt.show()

plt.plot(X, magsV, 'r+')
plt.xlabel('$1/\cos(z)$')
plt.ylabel('$Mag_V$')
plt.show()

#here is how you fit a line to data:
np.polyfit(X, magsB, 1)
np.polyfit(X, magsV, 1)
names = ['160269-V-1.fits', '160269-B-1.fits', '144284-V-1.fits', 
         '144284-B-1.fits']
dt = [2., 2., 0.5, 0.5]
K = [0.1839, 0.2467, 0.1839, 0.2467]
D = [48, 48, 55, 55]
M = [5, 5, 14, 14]
F = [-2.5 * np.log10(flux(photomDir + 'HD' + name)[0]/t) - k/np.cos(rad(d, m)) \
                                     for (name, t, k, d, m)\
                                     in zip(names,
                                            dt,
                                            K,
                                            M,
                                            D)]
F = np.array(F)

RMagB = [5.840, 5.840, 4.538, 4.538]
RMagV = [5.232, 5.232, 4.005, 4.005]
FLUXES = [flux(photomDir + 'HD' + name)[0] for name in names]
PMagB = [MAGB(aduv, adub, t, d, m) ]

def prevmag(name, tdp, serie, d, m):
    bFlux = flux(photomDir + 'HD' + name + '-B-' + str(serie) + '.fits')
    vFlux = flux(photomDir + 'HD' + name + '-V-' + str(serie) + '.fits')
    
    return MAGB(vFlux, bFlux, tdp, d, m), MAGV(vFlux, bFlux, tdp, d, m)
array([-0.43322667,  0.15047465,  0.6429345 , -0.21623384, -0.13389058,
        0.05429808, -0.13389461,  0.05429886])