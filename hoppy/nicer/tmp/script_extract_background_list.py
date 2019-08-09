#!/usr/bin/env python

import os 
import glob
import astropy.io.fits as pyfits
import numpy as np

"""
NICER spectrum
Channel 	Energy (keV)
0 			0.00-0.01 
1			0.01-0.02

	   (keV)		Channel
Band1: 0.2-0.4		20-39
Band2: 0.4-0.8		40-79
Band3: 0.8-1.6		80-159
Band4: 1.6-3.0		160-299
Badn5: 3.0-6.0		300-599
Band6: 6.0-12.0		600-1199
Band7: 12.0-15.0	1200-1499
"""

ebands = [[0.2,0.4],[0.4,0.8],[0.8,1.6],[1.6,3.0],[3.0,6.0],[6.0,12.0],[12.0,15.0]] # keV

def get_count(fitsfile,emin_keV,emax_keV,flag_rate=False):
	chmin = int(emin_keV*100)
	chmax = int(emax_keV*100)	
	hdu = pyfits.open(fitsfile)
	if flag_rate:
		exposure = hdu["SPECTRUM"].header['EXPOSURE']	
		flag_channel = np.logical_and(hdu["SPECTRUM"].data["CHANNEL"]>=chmin,hdu["SPECTRUM"].data["CHANNEL"]<chmax)
		count = sum(hdu["SPECTRUM"].data[flag_channel]["RATE"] * exposure)
	else:
		flag_channel = np.logical_and(hdu["SPECTRUM"].data["CHANNEL"]>=chmin,hdu["SPECTRUM"].data["CHANNEL"]<chmax)
		count = sum(hdu["SPECTRUM"].data[flag_channel]["COUNTS"])
	return count

def background_rate(obsid,flog):
	indir_list = glob.glob('out/*/%s' % obsid)
	print(indir_list)

	if len(indir_list) != 1:
		print("error : %s" % obsid)
		return -1

	indir = indir_list[0]
	src_pi = '%s/ni%s_nibackgen3C50_tot.pi' % (indir,obsid)
	bgd_pi = '%s/ni%s_nibackgen3C50_bgd.pi' % (indir,obsid)

	if not os.path.exists(src_pi) or not os.path.exists(bgd_pi):
		dump = "%s," % obsid
		dump += "skip,"
		dump += "\n"
		flog.write(dump)	
		return -1

	src_pi_hdu = pyfits.open(src_pi)
	bgd_pi_hdu = pyfits.open(bgd_pi)	
	exposure_src = float(src_pi_hdu["SPECTRUM"].header['EXPOSURE'])
	exposure_bgd = float(bgd_pi_hdu["SPECTRUM"].header['EXPOSURE'])
	mjdobs = src_pi_hdu["SPECTRUM"].header['MJD-OBS']

	dump  = "%s," % obsid
	dump += "%s," % src_pi_hdu["SPECTRUM"].header['OBJECT']
	dump += "%.5f," % mjdobs
	dump += "%s," % src_pi_hdu["SPECTRUM"].header['DATE-OBS']
	dump += "%s," % src_pi_hdu["SPECTRUM"].header['DATE-END']
	src_counts = []
	bgd_counts = []
	for eband in ebands:
		emin_keV = eband[0]
		emax_keV = eband[1]		
		src_counts.append(get_count(src_pi,emin_keV,emax_keV))
		bgd_counts.append(get_count(bgd_pi,emin_keV,emax_keV,flag_rate=True))
	dump += "%.3f," % exposure_src
	for i in range(len(ebands)):
		dump += "%d," % src_counts[i]
	for i in range(len(ebands)):
		dump += "%.3e," % (float(src_counts[i])/exposure_src)
	dump += "%.3f," % exposure_bgd
	for i in range(len(ebands)):
		dump += "%d," % bgd_counts[i]
	for i in range(len(ebands)):
		dump += "%.3e," % (float(bgd_counts[i])/exposure_bgd)		
	dump += '\n'
	print(dump)

	flog.write(dump)	

flog = open('bkgd_rxte_x.lst','w')
#num = 0
for obsid_path in glob.glob('out/BKGD_RXTE_?/*'):
	obsid = os.path.basename(obsid_path)
	background_rate(obsid,flog)
	#num += 1
	#if num > 10:
	#	break 

flog.close()	
