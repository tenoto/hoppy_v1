#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 September 13'
__version__ = '1.00'

import os 
import sys 
import argparse
import hoppy.nicer.nievent as nievt

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='nibarytime.py',
		usage='nibarytime.py infits ra dec orbfile',
		description='A barycentric correction script with adding barycentric time (MJD).',
		epilog='',
		add_help=True)

	parser.add_argument('infits',metavar='infits',type=str,        
		help='Input fits file.')
	parser.add_argument('ra',metavar='ra',type=float,        
		help='Input source RA (J2000).')	
	parser.add_argument('dec',metavar='dec',type=float,        
		help='Input source DEC (J2000).')		
	parser.add_argument('orbfile',metavar='orbfile',type=str,        
		help='Input orbital file.')	
	parser.add_argument('--outfits',metavar='outfits',type=str,        
		help='Output fits file.',default=None)
	parser.add_argument('--refframe',metavar='refframe',type=str,        
		help='Reference frame (default="ICRS").',default="ICRS")	
	parser.add_argument('--ephem',metavar='ephem',type=str,        
		help='Ephem (default="JPLEPH.430").',default="JPLEPH.430")							
	args = parser.parse_args()	
	print(args)

	nievt = nievt.NicerEventFits(args.infits)
	nievt.barycentric_correction(args.ra,args.dec,args.orbfile,
		outfits=args.outfits,refframe=args.refframe,ephem=args.ephem)
	
