#!/usr/bin/env python

import os 
import glob
import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy.table import Column, MaskedColumn

import matplotlib.pyplot as plt

from iminuit import Minuit
from probfit import Chi2Regression, linear

TessTimeBin_sec = 120.0 # sec
TessTimeBin_day = TessTimeBin_sec / 24. / 60. / 60. 

MISSING_VALUE = -9999

class Hist1D(object):

	def __init__(self, edges):
		self.edges = edges
		self.hist, edges = np.histogram([], bins=self.edges)
		self.bins = (edges[:-1] + edges[1:]) / 2.

	def fill(self, arr):
		hist, edges = np.histogram(arr, bins=self.edges)
		self.hist += hist

	@property
	def data(self):
		return self.bins, self.hist

def get_count_rate(cnt_array,exp_array):
	rate_list = []
	error_list = []
	for i in range(len(cnt_array.data[1])):
		cnt = cnt_array.data[1][i]
		exp = exp_array.data[1][i]

		if exp > 0:
			rate = float(cnt) / float(exp)
			error = float(np.sqrt(cnt)) / float(exp)
			#print(cnt,exp,rate)
		else:
			rate = 0
			error = 0
		rate_list.append(rate)
		error_list.append(error)
	return np.array(rate_list), np.array(error_list)

class TessLightCurve():
	def __init__(self,fitsfile):
		self.fitsfile = fitsfile
		print(self.fitsfile)

		self.hdu = fits.open(self.fitsfile)
		self.basename = os.path.splitext(os.path.basename(self.fitsfile))[0]
		print(self.basename)

		self.time_mjd = self.get_mjd()

		self.lc_orig_table = self.hdu['LIGHTCURVE'].data
		self.lc_orig_cols = self.lc_orig_table.columns

		self.edges = self.time_mjd+TessTimeBin_day/2.0 
		self.edges = np.insert(self.edges,0,self.time_mjd[0]-TessTimeBin_day/2.0)

		self.lc_list = {}

	def get_mjd(self):
		"""
		TUNIT1  = 'BJD - 2457000, days' / column units: Barycenter corrected TESS Julian
		TESS : BJD = TIME + 2457000 days 

		# MJD = BJD - 2400 000.5
		# https://en.wikipedia.org/wiki/Julian_day		
		"""				
		return self.hdu['LIGHTCURVE'].data['TIME'] + self.hdu['LIGHTCURVE'].header['BJDREFI'] + self.hdu['LIGHTCURVE'].header['BJDREFF']  - 2400000.5

	def cadence2mjd(self,a=0.00138893,b=58226.94810026):
		return a * self.hdu['LIGHTCURVE'].data['CADENCENO'] + b 

	def append_nicer_gti(self,input_niobs_list):
		self.niobsid_list = []
		self.nigti_list = []
		self.nimask = []

		for mjd in self.time_mjd:
			#out_gtinum = np.nan
			out_gtinum = MISSING_VALUE
			out_niobs = MISSING_VALUE
			out_mask = True
			for niobs in input_niobs_list:
				out_gtinum = niobs.get_mjdnum(mjd)
				#if not np.isnan(out_gtinum):
				if out_gtinum != MISSING_VALUE:
					out_niobs = niobs.obsid
					out_mask = False
					break
			#if not np.isnan(out_gtinum):
			#if out_gtinum != MISSING_VALUE:
			#	print(mjd,out_niobs,out_gtinum)
			self.niobsid_list.append(out_niobs)
			self.nigti_list.append(out_gtinum)
			self.nimask.append(out_mask)

	def append_nicer_count_rate(self,input_niobs_list,emin_keV,emax_keV):
		print(emin_keV)
		print(emax_keV)

		name_cnt = 'cnt_%s_%skeV' % (emin_keV,emax_keV)
		name_exp = 'exp_%s_%skeV' % (emin_keV,emax_keV)		
		name_rate = 'cps_%s_%skeV' % (emin_keV,emax_keV)				
		name_error = 'err_%s_%skeV' % (emin_keV,emax_keV)						

		lc_hist_cnt = Hist1D(edges=self.edges)
		lc_hist_exp = Hist1D(edges=self.edges)

		for niobs in input_niobs_list:
			#print(niobs.obsid)
			mask_energy = np.logical_and(niobs.keV>=emin_keV,niobs.keV<=emax_keV)
			lc_hist_cnt.fill(niobs.time_mjd[mask_energy])
			niobs.make_exposure_map()
			lc_hist_exp.fill(niobs.exposure_map)

		lc_hist_rate, lc_hist_error = get_count_rate(lc_hist_cnt,lc_hist_exp)

		self.lc_list[name_cnt] = lc_hist_cnt.data[1]
		self.lc_list[name_exp] = lc_hist_exp.data[1]
		self.lc_list[name_rate] = lc_hist_rate
		self.lc_list[name_error] = lc_hist_error

		#plt.step(*lc_hist_cnt.data)
		#plt.step(*lc_hist_exp.data)
		#plt.step(lc_hist_cnt.data[0],rate)
		#plt.errorbar(self.lc_hist_cnt.data[0],rate,
		#	yerr=error,marker='',drawstyle='steps-mid')
		#plt.xlim(58600-0.2,58600.3)
		#plt.savefig('test.pdf')

	def writeto(self,overwrite=True):
		outfitsfile = "%s_nicer.fits" % self.basename

		#new_column_nigti = np.array(self.nigti_list,dtype=np.int)
		#col = MaskedColumn(data=new_column_nigti,name='NIGTI',format='I',mask=self.nimask)
		#print(col)
		#exit()

		print(self.lc_list)

		new_column_mjd = np.array(self.time_mjd)
		new_column_niobsid = np.array(self.niobsid_list)
		new_column_nigti = np.array(self.nigti_list,dtype=np.int)
		new_columns = fits.ColDefs([
			fits.Column(name='MJD',format='F',array=new_column_mjd),			
			fits.Column(name='NIOBSID',format='10A',array=new_column_niobsid),
			#MaskedColumn(data=new_column_nigti,name='NIGTI',format='I',mask=self.nimask)			
			fits.Column(name='NIGTI',format='I',array=new_column_nigti)			
			])
		for keyword in self.lc_list:
			print(keyword,self.lc_list[keyword])
			new_columns += fits.ColDefs([
				fits.Column(name=keyword,format='F',array=self.lc_list[keyword])
				])
		hdu_primary = fits.PrimaryHDU()
		hdu_newtable = fits.BinTableHDU.from_columns(self.lc_orig_cols+new_columns,name='TRBLIST')	
		hdulist = fits.HDUList([hdu_primary,hdu_newtable])
		hdulist.writeto(outfitsfile,overwrite=True)	

class NicerEventFile():
	def __init__(self,fitsfile):
		self.fitsfile = fitsfile
		print(self.fitsfile)

		self.hdu = fits.open(self.fitsfile)
		self.obsid = str(self.hdu['EVENTS'].header['OBS_ID'])
		self.gti = self.hdu['GTI'].data
		self.gti_start_mjd = self.get_mjd(self.gti['START'])
		self.gti_stop_mjd = self.get_mjd(self.gti['STOP'])		
		self.time_mjd = self.get_mjd_events()

		self.keV = self.hdu['EVENTS'].data['PI'] / 100.0

	def get_mjd(self,timeseries):
		"""
		# MJD(TT) = (MJDREFI+MJDREFF) + (TIMEZERO+TIME)/86400
		# MJD(UTC) = (MJDREFI) + (TIMEZERO+TIME+LEAPINIT=2)/86400
		"""
		return self.hdu['EVENTS'].header['MJDREFI'] + self.hdu['EVENTS'].header['MJDREFF'] + (timeseries + self.hdu['EVENTS'].header['TIMEZERO']) / 86400.0

	def get_mjd_events(self):
		return self.get_mjd(self.hdu['EVENTS'].data['TIME'] )

	def get_mjdnum(self,mjd):
		#gtinum = np.nan		
		gtinum = MISSING_VALUE
		for i in range(len(self.gti)):
			#if self.gti_start_mjd[i] <= mjd < self.gti_stop_mjd[i]:
			#	gtinum = i
			if (mjd + TessTimeBin_day / 2.0 >= self.gti_start_mjd[i]) and (mjd - TessTimeBin_day / 2.0 <= self.gti_stop_mjd[i]):
				gtinum = i 		
				break 		
		return gtinum
		# 2789

	def show_mjd(self):
		for i in range(len(self.gti)):
			print(i,self.gti_start_mjd[i],self.gti_stop_mjd[i],(self.gti_stop_mjd[i]-self.gti_start_mjd[i])*24.*60.*60)

	def make_exposure_map(self,timebin_sec=1.0):
		self.exposure_map = []
		for i in range(len(self.gti)):
			start = self.gti_start_mjd[i]
			stop = self.gti_stop_mjd[i]
			time = start
			while time < stop:
				self.exposure_map.append(time)
				time += timebin_sec / 24. / 60. / 60. 

NicerEventFileList = []
for obsid_path in sorted(glob.glob('out/all/26050?????')):
	obsid = os.path.basename(obsid_path) 
	clevt = '%s/bary/ni%s_0mpu7_cl_bary.evt' % (obsid_path,obsid)
	if os.path.exists(clevt):
		NicerEventFileList.append(NicerEventFile(clevt))

data_dir_tess = "/Users/enoto/Dropbox/01_enoto/research/nicer/analysis/proximacen/data/200718_Notsu_TESS"
tess1 = TessLightCurve("%s/tess2019112060037-s0011-0000000388857263-0143-s_lc.fits" % data_dir_tess)
tess1.append_nicer_gti(NicerEventFileList)
tess1.append_nicer_count_rate(NicerEventFileList,0.2,0.3)
tess1.append_nicer_count_rate(NicerEventFileList,0.3,1.0)
tess1.append_nicer_count_rate(NicerEventFileList,1.0,3.0)
tess1.append_nicer_count_rate(NicerEventFileList,3.0,10.0)
tess1.append_nicer_count_rate(NicerEventFileList,10.0,15.0)
tess1.writeto()

data_dir_tess = "/Users/enoto/Dropbox/01_enoto/research/nicer/analysis/proximacen/data/200718_Notsu_TESS"
tess2 = TessLightCurve("%s/tess2019140104343-s0012-0000000388857263-0144-s_lc.fits" % data_dir_tess)
tess2.append_nicer_gti(NicerEventFileList)
tess2.append_nicer_count_rate(NicerEventFileList,0.2,0.3)
tess2.append_nicer_count_rate(NicerEventFileList,0.3,1.0)
tess2.append_nicer_count_rate(NicerEventFileList,1.0,3.0)
tess2.append_nicer_count_rate(NicerEventFileList,3.0,10.0)
tess2.append_nicer_count_rate(NicerEventFileList,10.0,15.0)
tess2.writeto()

