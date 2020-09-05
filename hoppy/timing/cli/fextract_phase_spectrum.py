#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import argparse

from hoppy.timing.timing import EventFits

__author__  = 'Teruaki Enoto'
__date__    = '2020 September 5'
__version__ = '0.01'

def get_parser():
	"""
	Creates a new argument parser.
	"""

	usage = """
plot event fits file (PULSE_PHASE)
"""

	parser = argparse.ArgumentParser(
		prog='fextract_phase_spectrum.py',
		usage='fextract_phase_spectrum.py infits filter_expression',
		description=usage,
		epilog='',
		add_help=True)
	version = '%(prog)s ' + __version__

	parser.add_argument('infits',metavar='infits',type=str,        
		help='Input event fits file.')	
	parser.add_argument('outpha',metavar='outpha',type=str,        
		help='Output phafile.')		
	parser.add_argument('filter_expression',metavar='filter_expression',type=str, default=None,
		help='filter expression')	
	parser.add_argument('--exposure_fraction',metavar='exposure_fraction',type=float, default=1.0,
		help='exposure fraction (default:1.0)')		
	return parser

def fextract_phase_spectrum(args):
	evtfits = EventFits(args.infits)
	evtfits.extract_phase_spectrum(args.outpha,args.filter_expression,
		exposure_fraction=args.exposure_fraction)
	print("finished...\n")

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	fextract_phase_spectrum(args)

if __name__=="__main__":
	main()

