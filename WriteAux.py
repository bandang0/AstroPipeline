#!/usr/bin/env python3
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds

from general import *

tdp = 0.5 # Temps de pose en seconde
serie = 1

#darkDir = dark2Dir # Tmp

darkPath = darkDir + 'mdark-{:.1f}'.format(tdp) + '.fits'
darkFits = fits.open(darkPath)
dfDark = pds.DataFrame(darkFits[0].data.flatten())
for filtre in ['B', 'V']:
    transPath = flatDir + 'mtrans-' + filtre + '.fits'
    transFits = fits.open(transPath)
    dfTrans = pds.DataFrame(transFits[0].data.flatten())
    fileList = glob.glob(HD144284Dir + 'HD144284-' + filtre + '-' + str (serie) +  '_*.fits')
    fileList = [f for f in fileList if not ".aux.fits" in f]
    for name in fileList:
        print(filtre, name)
        hduList = fits.open(name)
        hdu = hduList [0]
        dfObject = pds.DataFrame(hdu.data.flatten())
        hdu.data = ((dfObject - dfDark) / dfTrans).values.reshape(hdu.shape)
        hduList.writeto(name.split('.')[0] + '.aux.fits',
            output_verify = "ignore")

        hduList.close()
        print("done")
    transFits.close()
darkFits.close()
