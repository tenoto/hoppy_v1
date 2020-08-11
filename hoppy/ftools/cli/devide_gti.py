#!/usr/bin/env python

import os
import sys
import argparse
from argparse import ArgumentParser
import astropy.io.fits as fits

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-11 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('devide_gti.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
Devide the input event fitsfile into individual good time intervals. 
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--inputfile', '-i', type=str, 
		help='input event file.')	
	return parser

def devide_gti(inputfile,outdir='gti'):
	print("\n %s" % (sys._getframe().f_code.co_name))

	hdu = fits.open(inputfile)	
	extnames = []
	for i in range(len(hdu)):
		try:
			extnames.append(hdu[i].header['EXTNAME'])
		except:
			pass

	if 'GTI' in extnames:
		gti_extname = 'GTI'
		filetype = 'evt'
	elif 'STDGTI' in extnames:
		gti_extname = 'STDGTI'
		filetype = 'evt'	
	elif 'PREFILTER' in extnames:
		filetype = 'mkf'

	cmd  = 'rm -rf %s;mkdir -p %s\n' % (outdir,outdir)
	print(cmd);os.system(cmd)

	basename = os.path.splitext(os.path.basename(inputfile))[0]

	number_of_gti = len(hdu[gti_extname].data)
	for gtinum in range(number_of_gti):
		gti_dir = '%s/%003d' % (outdir,gtinum)
		cmd  = 'mkdir -p %s;\n' % gti_dir
		cmd += 'extract_gti_fits.py %s %d;\n' % (inputfile,gtinum)
		gti_basename = '%s_gti%03d' % (basename,gtinum)
		cmd += 'mv %s.{evt,gti} %s\n' % (gti_basename,gti_dir)
		print(cmd);os.system(cmd)

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	devide_gti(args.inputfile)

if __name__=="__main__":
	main()
