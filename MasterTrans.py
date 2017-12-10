#!/us/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import os
import pandas as pds

from general import *

filtre = "V"

flatDir = flat2Dir # Tmp

# Récupée fichies petinents
fileList = glob.glob (flatDir + 'Flat-' + filtre + '_*' + '.fits')

# On cée le maste dak
if len (fileList) == 0:
    print ("Pas de flat existant pour ces paramètres...")
else:
    # On lit un pemie fichie qui nous sevia de base
    print ("Reading : " + str (fileList [0].split (sep)[-1]))
    hduList = fits.open (fileList [0])
    hdu = hduList [0]
    header = hdu.header

    # On écupèe le temps de pose
    tdp = header ['Exp (ms)'] / 1000.0

    mdarkPath = darkDir + 'mdark-' + "{:.1f}".format(tdp) + '.fits'

    if os.path.isfile(mdarkPath):
        # On ouve le fichie mdak
        hduListMdark = fits.open (mdarkPath)
        hduMdark = hduListMdark [0]

        dfDark = pds.DataFrame (hduMdark.data.flatten ())
        dfFlat = pds.DataFrame (hdu.data.flatten())
        if tdp * 1000 != header ['Exp (ms)']:
            print ("Attention le temps de pose indiqué ne correspond pas à celui du fichier courant.")

        for i in range (1, len (fileList)):
            hduListCurrent = fits.open (fileList [i])
            print ("Reading : " + str (fileList [i].split (sep)[-1]))
            if tdp * 1000 != hduListCurrent [0].header ['Exp (ms)']:
                print ("Attention le temps de pose indiqué ne correspond pas à celui du fichier.")
            dfFlat [i] = hduListCurrent [0].data.flatten ()
            hduListCurrent.close ()

        dfFlat = pds.DataFrame (dfFlat.median (axis = 1))

        hdu.data = ((dfFlat - dfDark) / ((dfFlat - dfDark).mean ())).values.reshape (hdu.shape)

        hduList.writeto (flatDir + 'mtrans-' + filtre  + '.fits', output_verify = "ignore")
    else:
        print ("Fichier " + mdarkPath.split (sep)[-1] + " introuvable !")

    hduList.close ()
