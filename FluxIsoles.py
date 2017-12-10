#!/usr/bin/env python3

import math
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import numpy as np
import pandas as pds
import PythonPhot as pp

from general import *

#gaussB, psfB, psfmagB = pp.getpsf (imageB, xstarSumBV [ipal_B],  ypos,nmag,sky,1,1,np.arange(len(xpos)),5,'output_psf.fits')

def flux(name):
    print(name)
    image = fits.getdata (name)
    
    # On fit l'histogramme du fond de ciel avec une gaussienne
    skymod, skysig, skyskw = pp.mmm.mmm (image, readnoise = 0, minsky = 200)
    
    # On regarde à 5 sigma au dessus du fond de ciel pour identifier les étoiles
    hmin = skysig * 5
    
    # Largeur à mi-hauteur totale d'une étoile
    fwhm = 15
    
    # On repère les étoiles dans les images
    xstar, ystar, flux, sharp, round = pp.find.find (image, hmin, fwhm, [-2.0,2.0], [0.0,2.0])
   
    gain = 1.0
    
    # On utilise la fonction aper pour récupérer les magnitudes et les valeurs du fond de ciel pour les étoiles repérées
    mag, magerr, flux, fluxerr, sky, skyerr, badflag, outstr = \
            pp.aper.aper (image / gain, xstar , ystar, phpadu = gain, apr = 100, zeropoint = 25.11, skyrad = [5 * fwhm, 5 * fwhm], badpix = [-12000, 60000], exact = True, setskyval = 0)

    return np.nanmax(flux), np.nanmax(fluxerr), len (xstar)

# On récupère les valeurs de flux des étoiles de référence
    
outFilename = dataDir + str (dateObs2) + sep + "photom" + sep + "EtoilesReferences.log"
outFile = open (outFilename, "w")
    
obj = "HD202444"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "B"
numeroSerie = 1
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

outFile.write ("------------------ " + obj + " " + filtre + " " + str (numeroSerie) + " ------------------ \n\n")

l = []
for file in fileList:
    ret = flux (file)
    l.append (ret[0])
f  = np.mean (l)
outFile.write ("Flux calculés : " + str (np.array (l)) + "\n")
outFile.write ("Calcul de flux pour " + obj + " en filtre " + filtre + " - série " + str (numeroSerie) + " : " + str (f) + "\n\n")

obj = "HD202444"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "V"
numeroSerie = 1
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

outFile.write ("------------------ " + obj + " " + filtre + " " + str (numeroSerie) + " ------------------ \n\n")

l = []
for file in fileList:
    ret = flux (file)
    l.append (ret[0])
f  = np.mean (l)
outFile.write ("Flux calculés : " + str (np.array (l)) + "\n")
outFile.write ("Calcul de flux pour " + obj + " en filtre " + filtre + " - série " + str (numeroSerie) + " : " + str (f) + "\n\n")

obj = "HD203280"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "B"
numeroSerie = 1
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

outFile.write ("------------------ " + obj + " " + filtre + " " + str (numeroSerie) + " ------------------ \n\n")

l = []
for file in fileList:
    ret = flux (file)
    l.append (ret[0])
f  = np.mean (l)
outFile.write ("Flux calculés : " + str (np.array (l)) + "\n")
outFile.write ("Calcul de flux pour " + obj + " en filtre " + filtre + " - série " + str (numeroSerie) + " : " + str (f) + "\n\n")

obj = "HD203280"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "V"
numeroSerie = 1
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

outFile.write ("------------------ " + obj + " " + filtre + " " + str (numeroSerie) + " ------------------ \n\n")

l = []
for file in fileList:
    ret = flux (file)
    l.append (ret[0])
f  = np.mean (l)
outFile.write ("Flux calculés : " + str (np.array (l)) + "\n")
outFile.write ("Calcul de flux pour " + obj + " en filtre " + filtre + " - série " + str (numeroSerie) + " : " + str (f) + "\n\n")

obj = "HD214680"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "B"
numeroSerie = 1
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

outFile.write ("------------------ " + obj + " " + filtre + " " + str (numeroSerie) + " ------------------ \n\n")

l = []
for file in fileList:
    ret = flux (file)
    l.append (ret[0])
f  = np.mean (l)
outFile.write ("Flux calculés : " + str (np.array (l)) + "\n")
outFile.write ("Calcul de flux pour " + obj + " en filtre " + filtre + " - série " + str (numeroSerie) + " : " + str (f) + "\n\n")

obj = "HD214680"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "V"
numeroSerie = 1
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

outFile.write ("------------------ " + obj + " " + filtre + " " + str (numeroSerie) + " ------------------ \n\n")

l = []
for file in fileList:
    ret = flux (file)
    l.append (ret[0])
f  = np.mean (l)
outFile.write ("Flux calculés : " + str (np.array (l)) + "\n")
outFile.write ("Calcul de flux pour " + obj + " en filtre " + filtre + " - série " + str (numeroSerie) + " : " + str (f) + "\n\n")

outFile.close ()
    
