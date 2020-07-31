#!/usr/bin/env python

__author__ = 'Teruaki Enoto'
__version__ = '1.00'
# v1.00 : original version

import os 
import sys 
import argparse
import hoppy.xspec.xspec as xspec

if __name__=="__main__":
	parser = argparse.ArgumentParser(
		prog='get_xspec_rate.py',
		usage='get_xspec_rate.py pha rmf arf emin emax',
		description='Extract count rate in the specified energy band (emin <= E < emax).',
		epilog='',
		add_help=True)

	parser.add_argument('pha',metavar='pha',type=str,help='Input pha file.')
	parser.add_argument('rmf',metavar='rmf',type=str,help='Input rmf file.')	
	parser.add_argument('arf',metavar='arf',type=str,help='Input arf file.')		
	parser.add_argument('emin',metavar='emin',type=float,help='Minimum energy (keV).')        
	parser.add_argument('emax',metavar='emax',type=float,help='Maximum energy (keV).')        	
	args = parser.parse_args()      
	print(args)

	nipha = xspec.XspecPha(args.pha,backgrnd=None,
		rmffile=args.rmf,arffile=args.arf,
		outdir='tmp_xspec')
	print(nipha.get_rate_and_error(args.emin,args.emax))
