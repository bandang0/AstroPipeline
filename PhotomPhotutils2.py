#!/usr/bin/env python3

import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import datetime
import pandas as pds
import numpy as np
import sys

# Pour éviter les recursion error
sys.setrecursionlimit (10000)

from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.psf.sandbox import DiscretePRF
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.fitting import SLSQPLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IterativelySubtractedPSFPhotometry
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from general import *

# On lit les fichiers B et V
imageB, imageV = fits.getdata (photomDir + 'M13-B-1.aux.fits'), fits.getdata(photomDir + 'M13-V-1.aux.fits')
imageSumBV = fits.getdata (photomDir + 'M13-sumBV-1.aux.fits')

# Largeur à mi-hauteur totale d'une étoile
fwhm = 8
fshape = 17

# On commence par faire de la photométrie de base sur le fichier somme BV.
# Les données calculées pour le fichier somme serviront de point de départ à l'algorithme
# pour les images en filtre B et V.

################################## SUM BV ################################

# On commence par déterminer la déviation standard du fond de ciel.
bkgrms = MADStdBackgroundRMS ()
std = bkgrms (imageSumBV)

# On donne les paramètres du finder pour déterminer les étoiles de base dans l'image
iraffind = IRAFStarFinder (threshold = 10 * std,
                          fwhm = fwhm,
                          minsep_fwhm = 0.01, roundhi = 1.0, roundlo = -1.0,
                          sharplo = 0.1, sharphi = 0.8)



# On donne un critère de groupage des étoiles
daogroup = DAOGroup (2.0 * fwhm)

# On détermine le fond de ciel et la procédure de fitting
mmm_bkg = MMMBackground ()
#fitter = LevMarLSQFitter ()
fitter = SLSQPLSQFitter ()
psf_model = IntegratedGaussianPRF (sigma = fwhm / 2.35)

# On met tout ça dans une boîte noire qui fait des itérations soustractives
photometry = IterativelySubtractedPSFPhotometry (finder = iraffind,
                                                group_maker = daogroup,
                                                bkg_estimator = mmm_bkg,
                                                psf_model = psf_model,
                                                fitter = LevMarLSQFitter (),
                                                niters = 1, fitshape = (fshape, fshape))

# On exécute le tout et on extrait des résultats !
result_tabBV = photometry (image = imageSumBV)
residual_imageBV = photometry.get_residual_image ()
#result_tabBV.write ('result_tabBV.2.5s.dat', format = 'ascii')
#np.array (residual_imageBV).tofile ('residual_imageBV.2.5s.dat')
#residual_imageBV = np.fromfile ('residual_imageBV.2.5s.dat', dtype=np.float32).reshape ((1804, 1804))
#result_tabBV = astropy.table.table.Table ().read ('result_tabBV.2.5s.dat', format = 'ascii')

df = result_tabBV.to_pandas ()
df = df [(df ['x_0'] - 1150) ** 2 + (df ['y_0'] - 800) ** 2 > 600 ** 2]

result_tabBV = result_tabBV.from_pandas (df)

psf_modelBV = DiscretePRF.create_from_image (imageSumBV - MMMBackground ().calc_background (imageSumBV), result_tabBV, fshape, mode='median')
psf_modelB = DiscretePRF.create_from_image (imageB - MMMBackground ().calc_background (imageB), result_tabBV, fshape, mode='median')
psf_modelV = DiscretePRF.create_from_image (imageV - MMMBackground ().calc_background (imageV), result_tabBV, fshape, mode='median')

################################## SUM BV ################################

# On commence par déterminer la déviation standard du fond de ciel.
bkgrms = MADStdBackgroundRMS ()
std = bkgrms (imageSumBV)

# On donne les paramètres du finder pour déterminer les étoiles de base dans l'image
iraffind = IRAFStarFinder (threshold = 15 * std,
                          fwhm = fwhm,
                          minsep_fwhm = 0.01, roundhi = 1.0, roundlo = -1.0,
                          sharplo = 0.1, sharphi = 0.8)



# On donne un critère de groupage des étoiles
daogroup = DAOGroup (2.0 * fwhm)

# On détermine le fond de ciel et la procédure de fitting
mmm_bkg = MMMBackground ()
#fitter = LevMarLSQFitter ()
fitter = SLSQPLSQFitter ()
#psf_model = IntegratedGaussianPRF (sigma = fwhm / 2.35)

# On met tout ça dans une boîte noire qui fait des itérations soustractives
photometry = IterativelySubtractedPSFPhotometry (finder = iraffind,
                                                group_maker = daogroup,
                                                bkg_estimator = mmm_bkg,
                                                psf_model = psf_modelBV,
                                                fitter = LevMarLSQFitter (),
                                                niters = 1, fitshape = (fshape, fshape),
                                                aperture_radius = fwhm)

# On exécute le tout et on extrait des résultats !
result_tabBV = photometry (image = imageSumBV)
residual_imageBV = photometry.get_residual_image ()
result_tabBV.write ('result_tabBV.2.5s.dat', format = 'ascii')
np.array (residual_imageBV).tofile ('residual_imageBV.2.5s.dat')

# On ne garde que les colonnes utiles pour la suite (si on ne les supprime pas, ça fera planter les fonctions d'après)
result_tabBV.remove_column ('x_fit')
result_tabBV.remove_column ('y_fit')
result_tabBV.remove_column ('flux_0')
result_tabBV.remove_column ('flux_fit')
result_tabBV.remove_column ('id')
result_tabBV.remove_column ('group_id')
#result_tabBV.remove_column ('flux_unc')
#result_tabBV.remove_column ('x_0_unc')
#result_tabBV.remove_column ('y_0_unc')
result_tabBV.remove_column ('iter_detected')

#l1 = []
#l2 = []
#for i in range (2*fwhm):
#    for j in range (2*fwhm):
#        l1.append (i)
#        l2.append (j)
#r = psf_model.evaluate (np.array (l1), np.array (l2), 1.0, fwhm, fwhm)
#r = r.reshape ((2*fwhp, 2*fwhp))
#plt.imshow (r)

################################# FILTER B ################################

bkgrms = MADStdBackgroundRMS ()
std = bkgrms (imageB)

# Finder sans importance ici (sera ignoré) vu que l'on fournit les positions initiales détectées dans sumBV
iraffind = IRAFStarFinder (threshold = 2.5 * std,
                          fwhm = fwhm,
                          minsep_fwhm = 0.01, roundhi = 1.0, roundlo = -1.0,
                          sharplo = 0.1, sharphi = 0.8)

daogroup = DAOGroup (2.0 * fwhm)
mmm_bkg = MMMBackground ()
#fitter = LevMarLSQFitter ()
fitter = SLSQPLSQFitter ()
#psf_model = IntegratedGaussianPRF (sigma = fwhm / 2.35)

photometry = IterativelySubtractedPSFPhotometry (finder = iraffind,
                                                group_maker = daogroup,
                                                bkg_estimator = mmm_bkg,
                                                psf_model = psf_modelB,
                                                fitter = LevMarLSQFitter (),
                                                niters = 1, fitshape = (fshape, fshape),
                                                aperture_radius = fwhm)
#result_tabB = photometry (image = imageB, init_guesses = result_tabBV)
result_tabB = photometry (image = imageB)
residual_imageB = photometry.get_residual_image ()
result_tabB.write ('result_tabB.2.5s.dat', format = 'ascii')
np.array (residual_imageB).tofile ('residual_imageB.2.5s.dat')

################################# FILTER V ################################

bkgrms = MADStdBackgroundRMS ()
std = bkgrms (imageV)

# Finder sans importance ici (sera ignoré) vu que l'on fournit les positions initiales détectées dans sumBV
iraffind = IRAFStarFinder (threshold = 2.5 * std,
                          fwhm = fwhm,
                          minsep_fwhm = 0.01, roundhi = 1.0, roundlo = -1.0,
                          sharplo = 0.1, sharphi = 0.8)

daogroup = DAOGroup (2.0 * fwhm)
mmm_bkg = MMMBackground ()
#fitter = LevMarLSQFitter ()
fitter = SLSQPLSQFitter ()
#psf_model = IntegratedGaussianPRF (sigma = fwhm / 2.35)

photometry = IterativelySubtractedPSFPhotometry (finder = iraffind,
                                                group_maker = daogroup,
                                                bkg_estimator = mmm_bkg,
                                                psf_model = psf_modelV,
                                                fitter = LevMarLSQFitter (),
                                                niters = 1, fitshape = (fshape, fshape),
                                                aperture_radius = fwhm)
#result_tabV = photometry (image = imageV, init_guesses = result_tabBV)
result_tabV = photometry (image = imageV)
residual_imageV = photometry.get_residual_image ()
result_tabV.write ('result_tabV.2.5s.dat', format = 'ascii')
np.array (residual_imageV).tofile ('residual_imageV.2.5s.dat')

################################# PLOT #####################################

#norm = ImageNormalize (stretch = SqrtStretch ())
#plt.imshow (imageSumBV, cmap = 'Greys_r', origin = 'lower', norm = norm)
#plt.colorbar ()
#plt.plot (result_tabBV ['x_0'], result_tabBV ['y_0'], ls = 'none', color = 'cyan', marker = '+', ms = 10, lw = 1.5)

#plt.figure ()
#
#norm = ImageNormalize (stretch = SqrtStretch ())
#plt.imshow (imageB, cmap = 'Greys_r', origin = 'lower', norm = norm)
#plt.colorbar ()
#plt.plot (result_tabB [np.where (result_tabB ['iter_detected'] == 1)] ['x_0'], result_tabB [np.where (result_tabB ['iter_detected'] == 1)] ['y_0'], ls = 'none', color = 'cyan', marker = '+', ms = 10, lw = 1.5)
#
#plt.figure ()
#
#x = residual_imageB - np.min (residual_imageB)
##y = np.ma.masked_where(residual_imageB < 0, residual_imageB)
#y = residual_imageB
#plt.imshow (y, cmap = 'jet', vmax = 100., vmin=-100.)
#plt.colorbar ()

#
#norm = ImageNormalize (stretch = SqrtStretch ())
#plt.imshow (imageSumBV, cmap = 'Greys_r', origin = 'lower', norm = norm)
#plt.colorbar ()
#plt.plot (result_tabBV [np.where (result_tabBV ['iter_detected'] == 1)] ['x_0'], result_tabBV [np.where (result_tabBV ['iter_detected'] == 1)] ['y_0'], ls = 'none', color = 'cyan', marker = '+', ms = 10, lw = 1.5)
#plt.plot (result_tabBV [np.where (result_tabBV ['iter_detected'] == 2)] ['x_0'], result_tabBV [np.where (result_tabBV ['iter_detected'] == 2)] ['y_0'], ls = 'none', color = 'red', marker = '+', ms = 10, lw = 1.5)
#plt.plot (result_tabBV [np.where (result_tabBV ['iter_detected'] == 3)] ['x_0'], result_tabBV [np.where (result_tabBV ['iter_detected'] == 3)] ['y_0'], ls = 'none', color = 'green', marker = '+', ms = 10, lw = 1.5)
#plt.xlim (0, imageSumBV.shape [1] - 1)
#plt.ylim (0, imageSumBV.shape [0] - 1)
#
#x = np.sqrt (residual_imageBV - np.min (residual_imageBV))
#y = np.ma.masked_where(residual_imageBV < 0, residual_imageBV)
#plt.imshow ((y - y.min ())/(y.max () - y.min ()), cmap = 'jet', origin = 'lower')
#plt.colorbar ()
#plt.plot (result_tabBV [np.where (result_tabBV ['iter_detected'] == 1)] ['x_0'], result_tabBV [np.where (result_tabBV ['iter_detected'] == 1)] ['y_0'], ls = 'none', color = 'cyan', marker = '+', ms = 10, lw = 1.5)
#plt.plot (result_tabBV [np.where (result_tabBV ['iter_detected'] == 2)] ['x_0'], result_tabBV [np.where (result_tabBV ['iter_detected'] == 2)] ['y_0'], ls = 'none', color = 'red', marker = '+', ms = 10, lw = 1.5)
#plt.plot (result_tabBV [np.where (result_tabBV ['iter_detected'] == 3)] ['x_0'], result_tabBV [np.where (result_tabBV ['iter_detected'] == 3)] ['y_0'], ls = 'none', color = 'green', marker = '+', ms = 10, lw = 1.5)
#
