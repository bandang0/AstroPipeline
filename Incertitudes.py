#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:50:43 2017

@author: bandang0
"""

import math
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds
import PythonPhot as pp
import numpy as np

from general import *

#determine les incertitudes sur les flux mesurés des étoiles de reference 


starsB = [photomDir + i for i in ["HD160269-B-1.fits",
                                  "HD144284-B-1.fits",
                                  "HD147394-B-1.fits",
                                  "HD147394-B-2.fits",
                                  "HD147394-B-3.fits"]]

starsV = [photomDir + i for i in ["HD160269-V-1.fits",
                                  "HD144284-V-1.fits",
                                  "HD147394-V-1.fits",
                                  "HD147394-V-2.fits",
                                  "HD147394-V-3.fits"]]
serie = 1
dt = 0.5
filtre = 'B'
fileList = glob.glob(HD144284Dir + 'HD144284-' + filtre + '-' + str (serie) +  '_*.aux.fits')
fluxes = [flux(name)[0]/dt for name in fileList]

fluxes = np.array([float(i) for i in fluxes])
print(fluxes.std()/np.sqrt(10))
print(fluxes.std()/max(fluxes))                                  
