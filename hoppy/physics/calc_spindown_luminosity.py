#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 September 7'
__version__ = '1.00'

import argparse

import hoppy_physics 

if __name__=="__main__":

	description = 'A script to calculate spin-down luminosity from Period (s) and Pdot (s/s).'

	parser = argparse.ArgumentParser(
		prog='calc_spindown_luminosity.py',
		usage='calc_spindown_luminosity.py period pdot',
		description=description,
		epilog='',
		add_help=True,
		)

	parser.add_argument('period',metavar='period',type=float,help='Period (sec)')
	parser.add_argument('pdot',metavar='pdot',type=float,help='Period derivative (s/s)')	
	parser.add_argument('--moment_of_inertia',metavar='moment_of_inertia',type=float,
		default=1.0,help='moment_of_inertia (1e+45 g cm2) [default=1.0]')		
	args = parser.parse_args()	
	
	dump  = "Period: %.3e (s)\n" % args.period
	dump += "Pdot : %.3f (s/s)\n" % args.pdot
	dump += "moment_of_inertia: %.3e (g cm2)\n" % (args.moment_of_inertia * 1e+45)
	dump += "Spin-down luminosity: %.3e (erg/s)" % hoppy_physics.get_spindown_luminosity(args.period,args.pdot,moment_of_inertia=args.moment_of_inertia)
	print(dump)