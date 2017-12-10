#!/usr/bin/env python3

import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds
import numpy as np
import sys

from photutils.detection import IRAFStarFinder

from general import *

# On récupère les valeurs calculées par le script PhotomPhotutils
#result_tabBV = astropy.table.table.Table ().read ('result_tabBV.2.5s.dat', format = 'ascii')
result_tabV = astropy.table.table.Table ().read ('result_tabV.2.5s.dat', format = 'ascii')
result_tabB = astropy.table.table.Table ().read ('result_tabB.2.5s.dat', format = 'ascii')

df = result_tabV.to_pandas ()
df = df [(df ['x_fit'] - 1150) ** 2 + (df ['y_fit'] - 800) ** 2 > 400 ** 2]
result_tabV = result_tabV.from_pandas (df)

df = result_tabB.to_pandas ()
df = df [(df ['x_fit'] - 1150) ** 2 + (df ['y_fit'] - 800) ** 2 > 400 ** 2]
result_tabB = result_tabB.from_pandas (df)

# On fait attention au fait que les deux listes se réfèrent bien aux mêmes points
l = []
index1 = 0
index2 = 0
for elt1 in result_tabV:
    dMin = -1.
    indexSave = 0
    xCurr = elt1 ['x_fit']
    yCurr = elt1 ['y_fit']
    index2 = 0
    for elt2 in result_tabB:
        xCurr2 = elt2 ['x_fit']
        yCurr2 = elt2 ['y_fit']
        d = np.sqrt ((xCurr2 - xCurr) ** 2 + (yCurr2 - yCurr) ** 2)
        if d < dMin or dMin < 0:
            dMin = d
            indexSave = index2
        index2 += 1
    l.append ([index1, indexSave, dMin])
    index1 += 1
tab = pds.DataFrame (l)

# On garde que les points qui matchent bien à 5 pixels près
tab = tab [tab [2] < 5]
    
# On récupère les flux
fluxV = np.array(result_tabV ['flux_fit'].tolist ())
fluxB = np.array(result_tabB ['flux_fit'].tolist ())
#fluxVunc = np.array(result_tabV ['flux_unc'].tolist ())
#fluxBunc = np.array(result_tabB ['flux_unc'].tolist ())
#d = 

fluxB = fluxB [tab [1]]
fluxV = fluxV [tab [0]]
#fluxVunc = fluxVunc [tab [0]]
#fluxBunc = fluxBunc [tab [1]]
#
#relB = fluxBunc / fluxB
#relV = fluxVunc / fluxV
#
#fluxB = fluxB [(relB > 0.00001) & (relB < 1.) & (relV > 0.00001) & (relV < 1.)]
#fluxV = fluxV [(relB > 0.00001) & (relB < 1.) & (relV > 0.00001) & (relV < 1.)]

# On en tire les magnitudes
magV = MAGV (fluxV, fluxB, 30.0, 50, 59)
magB = MAGB (fluxV, fluxB, 30.0, 53, 46)

# On plot le diagramme
plt.plot (magB - magV, magV,'.',color='k')
plt.ylabel ('$M_V$')
plt.xlabel ('$M_B - M_V$')
plt.ylim ([18,11])
#plt.xlim ([0.3,0.6])
plt.figure ()
residual_imageB = np.fromfile ('residual_imageB.2.5s.dat', dtype=np.float32).reshape ((1804, 1804))
residual_imageV = np.fromfile ('residual_imageV.2.5s.dat', dtype=np.float32).reshape ((1804, 1804))

x = residual_imageV - np.min (residual_imageV)
#y = np.ma.masked_where(residual_imageB < 0, residual_imageB)
y = residual_imageB
plt.imshow (y, cmap = 'jet', vmax = 150., vmin=-150.)
plt.colorbar ()

