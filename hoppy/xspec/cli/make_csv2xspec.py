#!/usr/bin/env python

import argparse

import hoppy.xspec.xspec as xspec

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-12 : refactoring from a previous version  

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser(
		prog="make_csv2xspec.py",
		usage='%(prog)s filelist yamlfile outcsvfile',
		description="""
prepare csvfile for xspec fittings
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument(
		'filelist',metavar='filelist',type=str,
		help='observation file list.') 
	parser.add_argument(
		'yamlfile',metavar='yamlfile',type=str,
		help='parameter file in the yaml format.') 	
	parser.add_argument(
		'outcsvfile',metavar='outcsvfile',type=str,
		help='output csvfile') 		
	return parser

def make_csv2xspec(args):
	csv = xspec.CSVtoXSPEC(args.filelist,args.yamlfile)
	csv.make_csv2xspec(args.outcsvfile)

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	make_csv2xspec(args)

if __name__=="__main__":
	main()