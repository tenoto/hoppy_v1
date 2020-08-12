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
		prog="xspec_multi_fit.py",
		usage='%(prog)s  csvfile yamlfile',
		description="""
xspec fitting of multiple observations.
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument(
		'csvfile',metavar='csvfile',type=str,
		help='observation list in the csv format.') 
	parser.add_argument(
		'yamlfile',metavar='yamlfile',type=str,
		help='parameter file in the yaml format.') 	
	return parser

def xspec_multi_fit(args):
	manager = xspec.MonitoringManager(args.csvfile,args.yamlfile)
	manager.run_fit() 
	manager.collect_fitresults()

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	xspec_multi_fit(args)

if __name__=="__main__":
	main()