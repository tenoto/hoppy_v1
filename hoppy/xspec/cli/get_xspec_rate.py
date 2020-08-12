#!/usr/bin/env python

import os 
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
		prog="get_xspec_rate.py",
		usage='%(prog)s phafile emin emax [-b backgrnd] [-r rmffile] [-a arffile] ',
		description="""
Automatic xspec fitting.
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument(
		'phafile',metavar='phafile',type=str,
		help='source pha file for fitting.') 
	parser.add_argument(
		'emin',metavar='emin',type=float,
		help='Minimum energy (keV).')        
	parser.add_argument(
		'emax',metavar='emax',type=float,
		help='Maximum energy (keV).')   	
	parser.add_argument(
		'-b','--backgrnd',metavar='backgrnd',type=str,default=None,
		help='background pha file for fitting.') 		
	parser.add_argument(
		'-r','--rmffile',metavar='rmffile',type=str,default=None,
		help='rmffile for fitting.') 
	parser.add_argument(
		'-a','--arffile',metavar='arffile',type=str,default=None,
		help='arffile for fitting.') 
	return parser

def get_xspec_rate(args):
	cmd = 'rm -rf tmp_xspec_rate;\n'
	os.system(cmd)

	xspec_pha = xspec.XspecPha(args.phafile,
		outdir='tmp_xspec_rate',backgrnd=args.backgrnd,rmffile=args.rmffile,arffile=args.arffile)
	rate, rate_error = xspec_pha.get_rate_and_error(args.emin,args.emax)
	print("Rate: %.5f +/- %.5f (%.2f-%.2f keV)" % (rate,rate_error,args.emin,args.emax))

	cmd = 'rm -rf tmp_xspec_rate;\n'
	os.system(cmd)

	return rate, rate_error 

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	get_xspec_rate(args)

if __name__=="__main__":
	main()