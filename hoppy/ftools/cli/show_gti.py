#!/usr/bin/env python

import argparse
import astropy.io.fits as pyfits

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='show_gti.py',
		usage='show_gti.py input.fits',
		description='Show GTI information',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'input_fits',metavar='input_fits',type=str,        
		help='Input fits file.')
	args = parser.parse_args()	

	hdu = pyfits.open(args.input_fits)
	extname_list = []

	for i in range(len(hdu)):
		try:
			extname_list.append(hdu[i].header['EXTNAME'])
		except:
			pass
			#print('skip the extension...')

	if 'GTI' in extname_list:
		gtiname = 'GTI'
		filetype = 'evt'
	elif 'STDGTI' in extname_list:
		gtiname = 'STDGTI'
		filetype = 'evt'	
	elif 'PREFILTER' in extname_list:
		filetype = 'mkf'

	total_exp = 0.0

	if filetype == 'evt':
		num_of_gti = len(hdu[gtiname].data)
		print("# GTI-num: START STOP Exposure(s)")
		for i in range(num_of_gti):
			gti_start = hdu[gtiname].data[i]['START']
			gti_stop  = hdu[gtiname].data[i]['STOP']
			gti_exp   = gti_stop - gti_start
			print("%03d:  %.3f %.3f %.3f" % (i,gti_start,gti_stop,gti_exp))
			total_exp += gti_exp
	elif filetype == 'mkf':
		total_exp = len(hdu['PREFILTER'].data)
	print("Total exposure: %.3f (s)" % total_exp)

