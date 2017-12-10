#!/usr/bin/env python3
import math
import datetime
import os
import pandas as pds
from astropy.io import fits
import numpy as np

# Paramètes
if os.name == 'posix':
    dataDir = "/Users/bandang0/Google Drive/data/"
    sep = '/'
else:
    dataDir = "C:\\Users\\Bruno\\Google Drive OBS\\data\\"
    sep = '\\'


dateObs = datetime.date (2017, 10, 13)
dateObs2 = datetime.date (2017, 11, 22)
darkDir = dataDir + str(dateObs) + sep + 'darks' + sep
flatDir = dataDir + str(dateObs) + sep + 'flats' + sep
M13Dir = dataDir + str(dateObs) + sep + 'M13' + sep
HD144284Dir = dataDir + str(dateObs) + sep + 'HD144284' + sep
HD147394Dir = dataDir + str(dateObs) + sep + 'HD147394' + sep
HD160269Dir = dataDir + str(dateObs) + sep + 'HD160269' + sep
M39Dir = dataDir + str(dateObs2) + sep + 'M39' + sep
HD202444Dir = dataDir + str(dateObs2) + sep + 'HD202444' + sep
HD203280Dir = dataDir + str(dateObs2) + sep + 'HD203280' + sep
HD214680Dir = dataDir + str(dateObs2) + sep + 'HD214680' + sep
dark2Dir = dataDir + str(dateObs2) + sep + 'darks' + sep
flat2Dir = dataDir + str(dateObs2) + sep + 'flats' + sep
photomDir = dataDir + str(dateObs) + sep + 'photom' + sep

def df(path):
    return pds.DataFrame(fits.open(path)[0].data.flatten())

def rad(deg, minute):
    return math.pi*(deg + minute/60.0)/180.0

def mag(f, dt):
    return -2.5*np.log10(f/dt)

#Voici les coefficients photometriques ajustés pour les données de M13
# pour le 13 oct 2017
    
KB = 0.2467
KV = 0.1839
ZB = 21.53
ZV = 26.47
CB = 0.1848
CV = -9.889

def MAGV(ADUV, ADUB, tdp, degree, minute):
    # retourne la magnitude V dans le systeme de Oja et al. etant donné des
    # flux en ADU sur image traitée, le temps de pose en s et l'angle zenital
    # en degree et minute
    
    return (1/(1 + CV - CB)) * (CV * (mag(ADUB, tdp) + ZB - KB/np.cos(rad(degree, minute))) + (1 - CB) * (mag(ADUV, tdp) + ZV - KV/np.cos(rad(degree, minute))))
   
def MAGB(ADUV, ADUB, tdp, degree, minute):
    # retourne la magnitude B dans le systeme de Oja et al. etant donné des
    # flux en ADU sur image traitée, le temps de pose en s et l'angle zenital
    # en degree et minute
    
    return (1/(1 + CV - CB)) * ((1 + CV) * (mag(ADUB, tdp) + ZB - KB/np.cos(rad(degree, minute))) - CB * (mag(ADUV, tdp) + ZV - KV/np.cos(rad(degree, minute))))
  