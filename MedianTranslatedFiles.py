#!/usr/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds

from general import *

obj = "M39"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "B"
numeroSerie = 2

# Récupère fichiers pertinents
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str (numeroSerie) + '_*' + '_t.aux.fits')

# On crée le fichier médian
if len (fileList) == 0:
    print ("Pas de fichier existant pour ces paramètres...")
else:
    print ("Reading : " + str(fileList [0].split (sep)[-1]))
    hduList = fits.open (fileList [0])
    hdu = hduList [0]
    header = hdu.header

    dfTmp = pds.DataFrame (hdu.data.flatten())

    for i in range (1, len (fileList)):
        hduListCurrent = fits.open (fileList [i])
        print ("Reading : " + str (fileList [i].split (sep)[-1]))
        dfTmp [i] = hduListCurrent [0].data.flatten ()
        hduListCurrent.close ()

    hdu.data = dfTmp.median (axis = 1).values.reshape (hdu.shape)

    hduList.writeto (objDir + obj + "-" + filtre + "-" + str (numeroSerie) + '.fits', output_verify = "ignore")

    hduList.close ()
