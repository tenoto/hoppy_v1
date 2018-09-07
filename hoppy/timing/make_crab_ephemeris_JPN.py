#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 Sep. 7'
__version__ = '1.00'

code_description = """
aaa
"""

import argparse 

import make_crab_ephemeris 

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='make_crab_ephemeris_JPN.py',
		usage='make_crab_ephemeris_JPN.py param.yaml',
		description=code_description,
		epilog='',
		add_help=True,
		)
	parser.add_argument(
		'fyaml',metavar='fyaml',type=str,        
		help='yamlfile')
	args = parser.parse_args()	
	
	ephem = make_crab_ephemeris.Ephemeris_CrabPulsar(reference="JPN",fyaml=args.fyaml)



