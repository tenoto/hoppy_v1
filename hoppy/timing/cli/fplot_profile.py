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
		prog='fplot_profile.py',
		usage='fplot_profile.py infits',
		description=usage,
		epilog='',
		add_help=True)
	version = '%(prog)s ' + __version__

	parser.add_argument('infits',metavar='infits',type=str,        
		help='Input event fits file.')	
	parser.add_argument('--outpdf',metavar='outpdf',type=str, default="profile.pdf",  
		help='outpdf')	
	parser.add_argument('--phase_colname',metavar='phase_colname',type=str, default="PULSE_PHASE",  
		help='column name of the pulse phase (default: PULSE_PHASE)')	
	parser.add_argument('--nphase',metavar='nphase',type=int, default=20,  
		help='phase bin size (default 20)')	
											
	return parser

def fplot_profile(args):
	evtfits = EventFits(args.infits)
	x1,x1e,y1,y1e,ny1,ny1e = evtfits.get_profile(nphase=args.nphase,
		phase_colname=args.phase_colname)
	evtfits.plot_two_profiles(x1,x1e,y1,y1e,args.outpdf)
	#evtfits.plot_profile(
	#	outfits=args.outfits,
	#	phase_colname=args.phase_colname,
	#	binsize=args.binsize)
	print("finished...\n")

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	fplot_profile(args)

if __name__=="__main__":
	main()

