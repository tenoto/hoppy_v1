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
		prog="xspec_fit.py",
		usage='%(prog)s phafile [-b backgrnd] [-r rmffile] [-a arffile] [-m modelxcm] [-s binminsig] [-n binmaxbin]',
		description="""
Automatic xspec fitting.
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument(
		'phafile',metavar='phafile',type=str,
		help='source pha file for fitting.') 
	parser.add_argument(
		'-o','--outdir',metavar='outdir',type=str,default='out',
		help='output directory.') 		
	parser.add_argument(
		'-b','--backgrnd',metavar='backgrnd',type=str,default=None,
		help='background pha file for fitting.') 		
	parser.add_argument(
		'-r','--rmffile',metavar='rmffile',type=str,default=None,
		help='rmffile for fitting.') 
	parser.add_argument(
		'-a','--arffile',metavar='arffile',type=str,default=None,
		help='arffile for fitting.') 
	parser.add_argument(
		'-m','--modelxcm',metavar='modelxcm',type=str,default=None,
		help='modelxcm for fitting.') 				
	parser.add_argument(
		'-s','--binminsig',metavar='binminsig',type=int,default=5,
		help='binning minimum significance.') 					
	parser.add_argument(
		'-n','--binmaxbin',metavar='binmaxbin',type=int,default=50,
		help='binning max bins.') 
	parser.add_argument(
		'--fitemin',metavar='fitemin',type=float,default=0.4,
		help='fitting energy min (keV).') 								
	parser.add_argument(
		'--fitemax',metavar='fitemax',type=float,default=10.0,
		help='fitting energy max (keV).') 	
	parser.add_argument(
		'--ratebands',metavar='ratebands',type=str,default=[[0.4,6.0],[1.0,10.0]],
		help='rate energy bands (list) example:0.8-6.0,2.0-10.0 .') 											
	parser.add_argument(
		'--fluxbands',metavar='fluxbands',type=str,default=[[0.4,6.0],[1.0,10.0]],
		help='flux energy bands (list) example:0.8-6.0,2.0-10.0 .') 	
	parser.add_argument(
		'--parerrnum',metavar='parerrnum',type=str,default=[1,2,5],
		help='parameter error number list.') 

	return parser

def xspec_fit(args):
	xspec_pha = xspec.XspecPha(args.phafile,
		outdir=args.outdir,
		backgrnd=args.backgrnd,rmffile=args.rmffile,arffile=args.arffile,modelxcm=args.modelxcm,
		binminsig=args.binminsig,binmaxbin=args.binmaxbin,fitemin=args.fitemin,fitemax=args.fitemax,
		ratebands=args.ratebands,fluxbands=args.fluxbands,parerrnum=args.parerrnum)
	xspec_pha.run()

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	xspec_fit(args)

if __name__=="__main__":
	main()