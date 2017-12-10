#!/usr/bin/env python3

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds

from general import *

# Paramètes
tdp = 2.0
tdp1 = 1.5
tdp2 = 5.0

# Fits
hduList1 = fits.open ("mdark-" + "{:.1f}".format(tdp1) + '.fits')
hdu1 = hduList1 [0]
hduList2 = fits.open ("mdark-" + "{:.1f}".format(tdp2) + '.fits')
hdu2 = hduList2 [0]

# Intepolation
hdu1.data = hdu1.data + ((hdu2.data - hdu1.data) / (tdp2 - tdp1)) * (tdp - tdp1)
hdu1.header ['EXP (MS)'] = tdp * 1000

# Enegistement
hduList1.writeto ("mdark-" + "{:.1f}".format(tdp) + '.fits', output_verify = "ignore")

hduList1.close ()
hduList2.close ()
