#!/usr/bin/env python

import sys 
#import pyfits
import astropy.io.fits as pyfits

if len(sys.argv) != 2:
	sys.stderr.write("usage: %s infits\n" % sys.argv[0])
	quit()
infits = sys.argv[1]

total_exp = 0.0
hdu = pyfits.open(infits)
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

if filetype == 'evt':
	for line in hdu[gtiname].data:
		exp = line['STOP'] - line['START']
		#print(line['START'],line['STOP'],exp)
		total_exp += exp
	print(total_exp)
elif filetype == 'mkf':
	total_exp = len(hdu['PREFILTER'].data)
	print(total_exp)

