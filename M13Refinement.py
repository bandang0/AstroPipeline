#!/usr/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds
import numpy as np

from general import *

# On cherche à produire des fichiers fits permettant de faire de la photométrie avec PythonPhot.
# On produit :
# - un fichier fits sommant B et V translatés et recadré de manière à de pas avoir les bords du fichier
# - les fichiers individuels B et V recadrés de la même manière que le fichier somme

obj = "M13"
numeroSerie = 1
objB = photomDir + obj + "-B-" + str (numeroSerie) + ".fits"
objV = photomDir + obj + "-V-" + str (numeroSerie) + ".fits"

######################### PRECUT DES IMAGES B ET V ###########################

hduListB = fits.open (objB)
hduB = hduListB [0]

hduListV = fits.open (objV)
hduV = hduListV [0]

hduB.data = hduB.data [:1948, :1948]
hduB.header ['NAXIS1'] = 1948
hduB.header ['NAXIS2'] = 1948
hduB.header ['ROIWIDTH'] = 1948
hduB.header ['ROIHEIGH'] = 1948

hduV.data = hduV.data [100::, :1948]
hduV.header ['NAXIS1'] = 1948
hduV.header ['NAXIS2'] = 1948
hduV.header ['ROIWIDTH'] = 1948
hduV.header ['ROIHEIGH'] = 1948

hduListB.writeto (objB.replace (".fits", ".tmp.fits"), output_verify = "ignore")
hduListV.writeto (objV.replace (".fits", ".tmp.fits"), output_verify = "ignore")

hduListB.close ()
hduListV.close ()

############################ TRANSLATION ###########################

# On ouvre le fichier fits de référence
fileRef = objV.replace (".fits", ".tmp.fits")
hduListRef = fits.open (fileRef)
hduRef = hduListRef [0]
n = hduRef.header ['NAXIS1']
dataRef = hduRef.data

# On en extrait les pixels les plus significatifs (là où il y a des étoiles...)
indexList = []
threshold = pds.DataFrame (dataRef.flatten ()).quantile (0.995) [0]
for i in range (int (n / 4), int ( 3 * n / 4)):
    for j in range (int (n / 4), int ( 3 * n / 4)):
        if dataRef [i, j] > threshold:
            indexList.append ((i,j))

# On translate le fichier B
fileCurrent = objB.replace (".fits", ".tmp.fits")
hduListCurrent = fits.open (fileCurrent)
hduCurrent = hduListCurrent [0]
dataCurrent = hduCurrent.data

maxSum = -1.0
jOpt = -1
kOpt = -1

print ("Translations pour : " + str (fileCurrent.split (sep)[-1]))
print (hduCurrent.header['TIME'])

# On cherche le vecteur de translation optimal (celui qui maximise la corrélation)
for j in range (- 150, -120):
    for k in range (- 150, -120):
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
hduListCurrent.writeto (fileCurrent.replace (".fits", "_t.fits"), output_verify = "ignore")

#print (hduCurrent.data.flatten ().min ())

hduListCurrent.close ()

    
hduListRef.close ()

######################### CUT DES IMAGES B ET V ###########################

hduListB = fits.open (objB.replace (".fits", ".tmp_t.fits"))
hduB = hduListB [0]

hduListV = fits.open (objV.replace (".fits", ".tmp.fits"))
hduV = hduListV [0]

hduB.data = hduB.data [144:1948, 144:1948]
hduB.header ['NAXIS1'] = 1804
hduB.header ['NAXIS2'] = 1804
hduB.header ['ROIWIDTH'] = 1804
hduB.header ['ROIHEIGH'] = 1804

hduV.data = hduV.data [144::, 144::]
hduV.header ['NAXIS1'] = 1804
hduV.header ['NAXIS2'] = 1804
hduV.header ['ROIWIDTH'] = 1804
hduV.header ['ROIHEIGH'] = 1804

hduListB.writeto (objB.replace (".fits", ".aux.fits"), output_verify = "ignore")
hduListV.writeto (objV.replace (".fits", ".aux.fits"), output_verify = "ignore")

hduListB.close ()
hduListV.close ()

######################### CALCUL DE LA SOMME B + V #############################

hduListB = fits.open (objB.replace (".fits", ".aux.fits"))
hduB = hduListB [0]

hduListV = fits.open (objV.replace (".fits", ".aux.fits"))
hduV = hduListV [0]

hduB.data = hduB.data + hduV.data

hduListB.writeto (photomDir + obj + "-sumBV-" + str (numeroSerie) + ".aux.fits", output_verify = "ignore")

hduListB.close ()
hduListV.close ()

