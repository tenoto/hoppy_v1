#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script add "PULSE_PHASE" column calculated from 
# nu, nudot, nuddot, and epoch. 
# 
# HISTORY
# 2020-09-03: modify the structure of argpase
# 2018-09-14: Modified to compatible to photonphase (PINT) (Teru Enoto, ver 2.00)
# 2017-08-26: Translate to python script (Teru Enoto, ver 1.00)
# 2015-02-07: The original perl script aeplsphase.pl is written by Teruaki Enoto.

import os 
import argparse

from hoppy.timing.timing import EventFits

__author__  = 'Teruaki Enoto'
__date__    = '2020 September 3'
__version__ = '2.01'

def get_parser():
	"""
	Creates a new argument parser.
	"""

	usage = """
This python script calcurates pulse phases (PULSE_PHASE) of individual X-ray photons, based on an ephemeris of rotation frequency nu (Hz), its derivative (Hz/sec), and the folding epoch in a unit of TIME column. For negative optional value, use an equal symbol, e.g., --nudot=-3.687105e-10. We employed follogin formula. 

	nu*(t-t0) + nudot*(t-t0)^2/2.0 + nu2dot*(t-t0)^3/6.0 + nu3dot*(t-t0)^4/24.0 + nu4dot*(t-t0)^5/120.0

where t is an event time (TIME column), t0 is an epoch for folding. 

(Note) For folding with period and its derivative when constant P and Pdot are given, phases when each X-ray photons comes are calcurated by integrtion of 1/(P + t*Pdot) from 0 to t. Here t is a difference in time of the arrival time and the epoch. This gives, 1/(Pdot){ ln(1 + Pdot/P *t) }. 

(Reference); 
* http://heasarc.gsfc.nasa.gov/docs/xte/abc/time.html
* http://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
* http://www.starlink.rl.ac.uk/star/docs/sun67.htx/node226.html
"""

	parser = argparse.ArgumentParser(
		prog='faddphase_nu.py',
		usage='faddphase_nu.py infits epoch nu [--outfits] [--nudot] [--nudot2] [--nudot3] [--nudot4] [--offset]',
		description=usage,
		epilog='',
		add_help=True)
	version = '%(prog)s ' + __version__

	parser.add_argument('infits',metavar='infits',type=str,        
		help='Input event fits file.')	
	parser.add_argument('epoch',metavar='epoch',type=float,        
		help='Folding epoch in a unit of TIME.')
	parser.add_argument('nu',metavar='nu',type=float,        
		help='Folding frequency nu (Hz)')
	parser.add_argument('--nudot',metavar='nudot',type=float, default=0.0,  
		help='Folding frequency derivative nudot (Hz/s)')	
	parser.add_argument('--nu2dot',metavar='nu2dot',type=float, default=0.0,  
		help='Folding frequency 2nd derivative nu2dot (Hz/s2)')	
	parser.add_argument('--nu3dot',metavar='nu3dot',type=float, default=0.0,  
		help='Folding frequency 3rd derivative nu3dot (Hz/s3)')	
	parser.add_argument('--nu4dot',metavar='nu4dot',type=float, default=0.0,  
		help='Folding frequency 4th derivative nu4dot (Hz/s4)')	
	parser.add_argument('--outfits',metavar='outfits',type=str,default=None,
		help='Output fits file.')		
	parser.add_argument('--offset',metavar='offset',type=str,default="0.0",
		help='Optional offset to phase.')	
	parser.add_argument('--flag_mjd',action='store_true',dest='flag_mjd',
		default=False,help='Flag MJD (BARY_TIME) calculation.')  													
	return parser

def faddphase_nu(args):
	evtfits = EventFits(args.infits)

	evtfits.faddphase_nu(
		args.epoch,
		args.nu,
		nudot=float(args.nudot),
		nu2dot=float(args.nu2dot),
		nu3dot=float(args.nu3dot),
		nu4dot=float(args.nu4dot),
		outfits=args.outfits,
		offset=float(args.offset),
		flag_mjd=args.flag_mjd)	

	print("finished...\n")

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	faddphase_nu(args)

if __name__=="__main__":
	main()

