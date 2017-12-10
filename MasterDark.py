#!/usr/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds

from general import *

tdp = 6.5 # Temps de pose en seconde

# Récupére darks pertinents

darkDir = dark2Dir # Tmp

fileList = glob.glob (darkDir + 'Dark-' + '{:.1f}'.format(tdp) + '_*' + '.fits')

# On cée le maste dak
if len (fileList) == 0:
    print ("Pas de dark existant pour ces paramètres...")
else:
    print ("Reading : " + str(fileList [0].split (sep)[-1]))
    hduList = fits.open (fileList [0])
    hdu = hduList [0]
    header = hdu.header

    dfTmp = pds.DataFrame (hdu.data.flatten())
    if tdp * 1000 != header ['Exp (ms)']:
        print ("Attention le temps de pose indiqué ne correspond pas à celui du fichier.")

    for i in range (1, len (fileList)):
        hduListCurrent = fits.open (fileList [i])
        print ("Reading : " + str (fileList [i].split (sep)[-1]))
        if tdp * 1000 != hduListCurrent [0].header ['Exp (ms)']:
            print ("Attention le temps de pose indiqué ne correspond pas à celui du fichier.")
        dfTmp [i] = hduListCurrent [0].data.flatten ()
        hduListCurrent.close ()

    hdu.data = dfTmp.median (axis = 1).values.reshape (hdu.shape)

    hduList.writeto (darkDir + 'mdark-' + "{:.1f}".format(tdp) + '.fits', output_verify = "ignore")

    hduList.close ()
