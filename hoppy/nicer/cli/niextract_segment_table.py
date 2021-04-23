#!/usr/bin/env python

import os
import sys
import glob 
import argparse
import yaml
from astropy.io import fits
import pandas as pd
import numpy as np 

__author__ = 'Teruaki Enoto'
__version__ = '0.02'
# v0.02 : 2020-08-17 
# v0.01 : 2020-08-08 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('niextract_segment_table.py',
#		usage='%(prog)s ',
		description="""
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--inputdir', type=str, default=None,
		help='input directory path')						
	parser.add_argument('--outbasename', type=str, default=None,
		help='output basename')											
	return parser

def niextract_segment_table(args):
	print("\n %s" % (sys._getframe().f_code.co_name))
	print(args)

	inputdir = args.inputdir

	segment_df = pd.DataFrame({
		'OBJECT':[],
		'OBSID':[],
		'Segment': [],
		'Block':[],
		'EXPOSURE':[],
		'MJDOBS':[],
		'NCOUNTS':[],
		'Rate1':[],		
		'Rate1_error':[],				
		'Rate1_emin':[],
		'Rate1_emax':[],
		'Flux1':[],		
		'Flux1_errmin':[],
		'Flux1_errmax':[],								
		'Flux1_emin':[],		
		'Flux1_emax':[],				
		'MJD_DATEOBS':[],
		'MJD_DATEEND':[],
		'DATEOBS':[],
		'DATEEND':[],
		'TSTART':[],
		'TSTOP':[],
		})

	target = '%s/proc/*/segment/seg*/block/block*/fit/ni*_0mpu7_seg*_3c50_tot_*_fit.yaml' % (inputdir)
	list_fit_yamlfile = glob.glob(target)
	for fit_yamlfile in list_fit_yamlfile:
		#print(fit_yamlfile)

		segment = fit_yamlfile.split('/segment/')[-1].split('/block/')[0]
		block = fit_yamlfile.split('/block/')[-1].split('/fit/')[0]

		param = yaml.load(open(fit_yamlfile),Loader=yaml.FullLoader)

		hdu = fits.open(param['phafile'])
		emin_keV = float(param['ratebands'][0][0])
		emax_keV = float(param['ratebands'][0][1])		

		flag = np.logical_and(
			emin_keV * 1000 < (hdu['SPECTRUM'].data['CHANNEL'] * 10),
			(hdu['SPECTRUM'].data['CHANNEL'] * 10) <= emax_keV * 1000
			)
		ncounts = sum(hdu['SPECTRUM'].data['COUNTS'][flag])

		segment_df = segment_df.sort_values(by=['MJDOBS'])
		segment_df = segment_df.append({
			'OBJECT': param['OBJECT'],			
			'OBSID': param['OBSID'],
			'Segment': segment,
			'Block': block,
			'EXPOSURE': param['EXPOSURE'],
			'MJDOBS': param['MJDOBS'],
			'NCOUNTS':ncounts, 
			'Rate1': param['ratelist'][0][0],
			'Rate1_error': param['ratelist'][0][1],
			'Rate1_emin': param['ratebands'][0][0],
			'Rate1_emax': param['ratebands'][0][1],
			'Flux1': param['fluxlist'][0][0],
			'Flux1_errmin': param['fluxlist'][0][1],
			'Flux1_errmax': param['fluxlist'][0][2],
			'Flux1_emin': param['fluxbands'][0][0],									
			'Flux1_emax': param['fluxbands'][0][1],												
			'MJD_DATEOBS': param['MJD_DATEOBS'],			
			'MJD_DATEEND': param['MJD_DATEEND'],			
			'DATEOBS': param['DATEOBS'],			
			'DATEEND': param['DATEEND'],			
			'TSTART': param['TSTART'],															
			'TSTOP': param['TSTOP'],																		
			}, 
			ignore_index=True)
#		if param['ratelist'][0][0] < 0.8:
#			print("***",param['OBSID'],segment,block,param['ratelist'][0][0])
	print(segment_df)

	outbasename = '%s/nipipeline_block' % args.inputdir
	csvfile = outbasename + '.csv'
	htmlfile = outbasename + '.html'	
	segment_df.to_csv(csvfile)
	segment_df.to_html(open(htmlfile,'w'))

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	niextract_segment_table(args)

if __name__=="__main__":
	main()