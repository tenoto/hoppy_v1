#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 September 7'
__version__ = '1.00'

import argparse

import hoppy_physics 

if __name__=="__main__":

	description = 'A script to calculate luminosity from flux and distance.'

	parser = argparse.ArgumentParser(
		prog='calc_luminosity.py',
		usage='calc_luminosity.py flux d_kpc',
		description=description,
		epilog='',
		add_help=True,
		)

	parser.add_argument('flux',metavar='flux',type=float,help='flux (erg/s/cm2)')
	parser.add_argument('d_kpc',metavar='d_kpc',type=float,help='distance (kpc)')	
	args = parser.parse_args()	
	
	dump  = "flux: %.3e (erg/s/cm2)\n" % args.flux
	dump += "distance : %.3f (kpc)\n" % args.d_kpc
	dump += "luminosity: %.3e (erg/s)" % get_luminosity(args.flux,args.d_kpc)
	print(dump)