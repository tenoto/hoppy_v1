#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 September 8'
__version__ = '1.00'

import argparse

import hoppy_physics 

if __name__=="__main__":

	description = 'A script to calculate survival function, 1-cdf, at a specified value.'

	parser = argparse.ArgumentParser(
		prog='calc_gaussian_sf.py',
		usage='calc_gaussian_sf.py mean std value',
		description=description,
		epilog='',
		add_help=True,
		)

	parser.add_argument('mean',metavar='mean',type=float,help='Mean of the Gaussian distribution.')
	parser.add_argument('std',metavar='std',type=float,help='Standard deviation of the distribution.')	
	parser.add_argument('value',metavar='value',type=float,help='Evaluated value of the distribution.')		
	args = parser.parse_args()	
	
	dump  = "mean: %.3e \n" % args.mean
	dump += "std : %.3e \n" % args.std
	dump += "valu: %.3e \n" % args.value
	dump += "sf: %.4e" % hoppy_physics.get_norm_sf(args.mean,args.std,args.value)
	print(dump)