#!/usr/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds
import numpy as np
import PythonPhot as pp

from general import *

# On lit les fichiers B et V
imageB, imageV = fits.getdata (photomDir + 'M13-B-1.aux.fits'), fits.getdata(photomDir + 'M13-V-1.aux.fits')
imageSumBV = fits.getdata (photomDir + 'M13-sumBV-1.aux.fits')

# On fit l'histogramme du fond de ciel avec une gaussienne
skymodB, skysigB, skyskwB = pp.mmm.mmm (imageB, readnoise = 0, minsky = 200)
skymodV, skysigV, skyskwV = pp.mmm.mmm (imageV, readnoise = 0, minsky = 200)
skymodSumBV, skysigSumBV, skyskwSumBV = pp.mmm.mmm (imageSumBV, readnoise = 0, minsky = 200)

# On regarde à 5 sigma au dessus du fond de ciel pour identifier les étoiles
hminB = skysigB * 2
hminV = skysigV * 2
hminSumBV = skysigSumBV * 2

# Largeur à mi-hauteur totale d'une étoile
fwhm = 12

# On repère les étoiles dans les images
#xstarB, ystarB, fluxB, sharpB, roundB = pp.find.find (imageB, hminB, fwhm, [-1.0,1.0], [0.2,1.0])
#xstarV, ystarV, fluxV, sharpV, roundV = pp.find.find (imageV, hminV, fwhm, [-1.0,1.0], [0.2,1.0])
xstarSumBV, ystarSumBV, fluxSumBV, sharpSumBV, roundSumBV = pp.find.find (imageSumBV, hminSumBV, fwhm, [-1.0,1.0], [0.1,0.8])

# On affiche ce qu'on a préparé
plt.ion ()
plt.imshow (imageSumBV, vmin = 0, vmax = 20 * skysigSumBV, cmap = 'Greys_r')
plt.plot ( xstarSumBV, ystarSumBV, 'o', ms = 5, mfc = 'none', lw = 2, mec = 'r')
#plt.xlim ([400,1350]); plt.ylim([450,1250])

gainB = 1.0
gainV = 1.0

# On ne considère qu'une portion de l'image pour bien isoler M13
ipal_B = ((xstarSumBV > 603) & (xstarSumBV < 2048) &
          (xstarSumBV > 0) & (xstarSumBV < 1350))
ipal_V = ((xstarSumBV > 603) & (xstarSumBV < 2048) &
          (xstarSumBV > 0) & (xstarSumBV < 1350))

# On utilise la fonction aper pour récupérer les magnitudes et les valeurs du fond de ciel pour les étoiles repérées
magB, magerrB, fluxB, fluxerrB, skyB, skyerrB, badflag, outstr = \
        pp.aper.aper (imageB / gainB, xstarSumBV [ipal_B], ystarSumBV [ipal_B], phpadu = gainB, apr = 5, zeropoint = 25.11, skyrad = [3 * fwhm, 5 * fwhm], badpix = [-12000, 60000], exact = True, setskyval = 0)
magV, magerrV, fluxV, fluxerrV, skyV, skyerrV, badflag, outstr = \
        pp.aper.aper (imageV / gainV, xstarSumBV [ipal_V], ystarSumBV [ipal_V], phpadu = gainV, apr = 5, zeropoint = 25.11, skyrad = [3 * fwhm, 5 * fwhm], badpix = [-12000, 60000], exact = True, setskyval = 0)

# On détermine la fonction de PSF globale des étoiles détectées (getpsf fait juste un fit sur la première des étoiles détectées)
gaussB, psfB, psfmagB = pp.getpsf.getpsf (imageB, xstarSumBV [ipal_B],  ystarSumBV [ipal_B], magB, skyB, 1, 1, np.arange (len (xstarSumBV [ipal_B])), fwhm + 1, fwhm, 'M13-B-1-getpsf.aux.fits')
gaussV, psfV, psfmagV = pp.getpsf.getpsf (imageV, xstarSumBV [ipal_V],  ystarSumBV [ipal_V], magV, skyV, 1, 1, np.arange (len (xstarSumBV [ipal_V])), fwhm + 1, fwhm, 'M13-V-1-getpsf.aux.fits')

# On groupe les étoiles selon leur distances relatives
ngroupB = pp.group.group (xstarSumBV [ipal_B],  ystarSumBV [ipal_B], fwhm)
ngroupV = pp.group.group (xstarSumBV [ipal_V],  ystarSumBV [ipal_V], fwhm)

# On fait un fit simultané de PSF


plt.plot (magB - magV, magV,'.',color='k')
plt.ylim ([15,13])
#plt.xlim ([0.3,0.6])
plt.ylabel ('$magV$')
plt.xlabel ('$magB - magV$')
