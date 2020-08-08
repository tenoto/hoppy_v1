#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from astropy.time import Time

TessTimeBin = 120.0 # sec



class TessLightCurve():
	def __init__(self,fitsfile):
		self.fitsfile = fitsfile
		print(self.fitsfile)

		self.hdu = fits.open(self.fitsfile)

		# TUNIT1  = 'BJD - 2457000, days' / column units: Barycenter corrected TESS Julian
		# TESS : BJD = TIME + 2457000 days 
		#			 = 

class NicerEventFile():
	def __init__(self,fitsfile):
		self.fitsfile = fitsfile
		print(self.fitsfile)

		self.hdu = fits.open(self.fitsfile)

		# MJD = BJD - 2400 000.5
		# BJD: TESS = NICER + 2400 000.5 - 2457 000.0 = NICER -56,999.5
		# MJD(TT) = (MJDREFI+MJDREFF) + (TIMEZERO+TIME)/86400
		# MJD(UTC) = (MJDREFI) + (TIMEZERO+TIME+LEAPINIT=2)/86400

data_dir_tess = "/Users/enoto/Dropbox/01_enoto/research/nicer/analysis/proximacen/data/200718_Notsu_TESS"
tess1 = TessLightCurve("%s/tess2019112060037-s0011-0000000388857263-0143-s_lc.fits" % data_dir_tess)
tess2 = TessLightCurve("%s/tess2019140104343-s0012-0000000388857263-0144-s_lc.fits" % data_dir_tess)
print(tess1.hdu['LIGHTCURVE'].data['TIME'])
print(tess2.hdu['LIGHTCURVE'].data['TIME'])
print(tess2.hdu['LIGHTCURVE'].data['CADENCENO'])

print(max(np.abs(tess2.hdu['LIGHTCURVE'].data['TIME']-(120.0/24/60/60*tess2.hdu['LIGHTCURVE'].data['CADENCENO'] + 1227.4594823256))))

"""
data_dir_nicer = "/Users/enoto/work/project/hoppy_nicer/out/proximacen/"
obsid = 2605012003 
nicer = NicerEventFile("%s/%s/bary/ni%s_0mpu7_cl_bary.evt" % (data_dir_nicer,obsid,obsid))
print(nicer.hdu['EVENTS'].data['TIME'])
print(nicer.hdu['EVENTS'].header['MJDREFI'])
print(nicer.hdu['EVENTS'].header['MJDREFF'])
print(nicer.hdu['EVENTS'].header['DATE-OBS'])

time_nicer2tess = nicer.hdu['EVENTS'].data['TIME']/86400.0 + nicer.hdu['EVENTS'].header['MJDREFI'] + nicer.hdu['EVENTS'].header['MJDREFF'] - 56999.5
print(time_nicer2tess)
"""