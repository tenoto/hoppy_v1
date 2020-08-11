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
	parser = argparse.ArgumentParser('niextract_spec.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
Extract spectrum of the NICER data.
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--inputfile', '-i', type=str, 
		help='input event file.')
	parser.add_argument('--outfile', '-o', type=str, 
		help='output pha file.')	
	parser.add_argument('--exclude_detid', '-e', type=str, default="14,34", 
		help='extclude det_id (default="14,34").')			
	return parser

def niextract_spec(inputfile,outfile,exclude_detid):
	print("...%s" % (sys._getframe().f_code.co_name))

	exclude_detid_list = exclude_detid.split(',')
	filter_expr = ''
	for i in range(len(exclude_detid_list)):
		if i > 0:
			filter_expr += " .AND. "
		filter_expr += "( DET_ID != %s )" % exclude_detid_list[i]
	print(filter_expr)		
	cmd = 'rm -f tmp_niextract_spec.evt\n'
	print(cmd);os.system(cmd)
	cmd = 'fselect %s tmp_niextract_spec.evt "%s"' % (inputfile,filter_expr)
	print(cmd);os.system(cmd)

	cmd  = 'xselect<<EOF\n'
	cmd += 'xsel\n'
	cmd += 'read event tmp_niextract_spec.evt ./\n' 
	cmd += 'yes\n'
	cmd += 'extract spec\n'
	cmd += 'save spec %s\n' % outfile
	cmd += 'yes\n'
	cmd += 'exit\n'
	cmd += 'no\n'
	cmd += 'exit\n'	
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)		

	cmd = 'rm -f tmp_niextract_spec.evt xselect.log\n'
	print(cmd);os.system(cmd)	

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	niextract_spec(args.inputfile,args.outfile,args.exclude_detid)

if __name__=="__main__":
	main()
