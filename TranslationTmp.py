#!/usr/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob
import datetime
import pandas as pds

from general import *

obj = "M13"
objDir = dataDir + str(dateObs) + sep + obj + sep
filtre = "V"
numeroSerie = 1

# Récupère fichiers pertinents
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

# On ouvre le fichier fits de référence
fileRef = fileList [0]
hduListRef = fits.open (fileRef)
hduRef = hduListRef [0]
n = hduRef.header ['NAXIS1']
dataRef = hduRef.data.flatten ()

# On en extrait les pixels les plus significatifs (là où il y a des étoiles...)
indexList = []
threshold = pds.DataFrame (dataRef).quantile (0.9999) [0]
for i in range (0, n * n):
    if dataRef [i] > threshold:
        indexList.append (i)

# On parcourt tous les fichiers
for i in range (1, len (fileList)):
    hduListCurrent = fits.open (fileList [i])
    hduCurrent = hduListCurrent [0]
    dataCurrent = hduCurrent.data.flatten ()
    
    maxSum = -1.0
    tOpt = 0
    jOpt = 0
    kOpt = 0
    
    print ("Translations pour : " + str (fileList [i].split (sep)[-1]))
                
    # On utilise la valeur de translation repérée et on enregistre le fichier résultant
    dataFinal = np.zeros (len (dataCurrent), dtype=np.float32)
    print ("     Translation optimale : j = " + str (jOpt) + " et k = " + str (kOpt))
    print ("     MaxSum = " + str (maxSum))
    for l in range (0, n * n):
        dataFinal [l] = abs (dataCurrent [l] - dataRef [l])
    hduCurrent.data = dataFinal.reshape (hduCurrent.shape)
    hduListCurrent.writeto (fileList [i].replace (".aux.fits", "_t.aux.tmp.fits"), output_verify = "ignore")
    
    #print (hduCurrent.data.flatten ().min ())
    
    hduListCurrent.close ()
    
    
hduListRef.close ()
