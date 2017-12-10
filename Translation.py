#!/usr/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob
import datetime
import pandas as pds

from general import *

obj = "M39"
objDir = dataDir + str(dateObs2) + sep + obj + sep
filtre = "V"
numeroSerie = 2

# Récupère fichiers pertinents
fileList = glob.glob (objDir + obj + "-" + filtre + "-" + str(numeroSerie) + '_*' + 'aux.fits')
fileList = [f for f in fileList if not "_t.aux.fits" in f]

# On ouvre le fichier fits de référence
fileRef = fileList [0]
hduListRef = fits.open (fileRef)
hduRef = hduListRef [0]
n = hduRef.header ['NAXIS1']
dataRef = hduRef.data

# On en extrait les pixels les plus significatifs (là où il y a des étoiles...)
indexList = []
threshold = pds.DataFrame (dataRef.flatten ()).quantile (0.999) [0]
for i in range (int (n / 4), int ( 3 * n / 4)):
    for j in range (int (n / 4), int ( 3 * n / 4)):
        if dataRef [i, j] > threshold:
            indexList.append ((i,j))

# On parcourt tous les fichiers
for i in range (1, len (fileList)):
    hduListCurrent = fits.open (fileList [i])
    hduCurrent = hduListCurrent [0]
    dataCurrent = hduCurrent.data
    
    maxSum = -1.0
    jOpt = -1
    kOpt = -1
    
    print ("Translations pour : " + str (fileList [i].split (sep)[-1]))
    print (hduCurrent.header['TIME'])
    
    # On cherche le vecteur de translation optimal (celui qui maximise la corrélation)
    for j in range (- int (n / 16), int (n / 16)):
        for k in range (- int (n / 16), int (n / 16)):
            sumTmp = 0.0
            for l in indexList: # Pour chaque valeur du vecteur de translation, on calcule la corrélation sur les pixels significatifs de l'image de référence 
                sumTmp = sumTmp + dataRef [l] * dataCurrent [l[0] + j, l[1] + k]
            if sumTmp > maxSum: # Si on trouve une corrélation plus élevée que précédemment, on enregistre le résultat
                maxSum = sumTmp
                jOpt = j
                kOpt = k
                
    # On utilise la valeur de translation repérée et on enregistre le fichier résultant
    dataFinal = np.zeros ((n, n), dtype=np.float32)
    print ("     Translation optimale : j = " + str (jOpt) + " et k = " + str (kOpt))
    print ("     MaxSum = " + str (maxSum))
    for j in range (0, n):
        for k in range (0, n):
            if (j + jOpt >= 0 and j + jOpt < n and k + kOpt >= 0 and k + kOpt < n):
                dataFinal [j, k] = dataCurrent [j + jOpt, k + kOpt]
    hduCurrent.data = dataFinal
    hduListCurrent.writeto (fileList [i].replace (".aux.fits", "_t.aux.fits"), output_verify = "ignore")
    
    #print (hduCurrent.data.flatten ().min ())
    
    hduListCurrent.close ()
    
    
hduListRef.close ()
