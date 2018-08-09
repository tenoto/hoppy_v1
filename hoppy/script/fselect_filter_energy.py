#!/usr/bin/env python

__author__  = 'Teruaki Enoto'
__version__ = '1.01'
__date__    = '2018 August 9'

import os 
import sys 
import argparse
import astropy.io.fits as pyfits

parser = argparse.ArgumentParser(
	prog='fselect_filter_energy.py',
	usage='fselect_filter_energy.py inputevtfits outputevtfits emin emax',
	description='filter an event fit file with energy.',
	epilog='',
	add_help=True,
	)

parser.add_argument(
	'inputevtfits',metavar='inputevtfits',type=str,        
	help='Input event fits filename.')
parser.add_argument(
	'outputevtfits',metavar='outputevtfits',type=str,        
	help='Output event fits filename.')
parser.add_argument(
	'emin',metavar='emin',type=float,        
	help='Energy min (keV) for the event filtering.')
parser.add_argument(
	'emax',metavar='emax',type=float,        
	help='Energy max (keV) for the event filtering.')	
args = parser.parse_args()	


sys.stdout.write('input event fits file : %s\n' % args.inputevtfits)	
sys.stdout.write('output event fits file : %s\n' % args.outputevtfits)	
sys.stdout.write('energy min (keV): %.2f\n' % args.emin)
sys.stdout.write('energy max (keV): %.2f\n' % args.emax)

if not os.path.exists(args.inputevtfits):
	sys.stderr.write('input event file does not exist. %s \n' % args.inputevtfits)
	exit()
if os.path.exists(args.outputevtfits):
	sys.stderr.write('output event file has already existed. %s \n' % args.outputevtfits)
	exit()	

hdu = pyfits.open(args.inputevtfits)
if len(hdu['EVENTS'].data) == 0:
	sys.stdout.write('Skip: No events in fits file %s' % args.inputevtfits)
	quit()
telescope  = hdu[0].header['TELESCOP']
instrument = hdu[0].header['INSTRUME']
KEV_TO_PI = None
if telescope == 'NICER' and instrument == 'XTI':
	KEV_TO_PI = 1000.0/10.0 # Conversion from PI to keV (PI is in units of 10 eV)
sys.stdout.write('TELESCOP: %s\n' % telescope)
sys.stdout.write('INSTRUME: %s\n' % instrument)
sys.stdout.write('KEV_TO_PI: %.3f\n' % KEV_TO_PI) 

pi_min = int(KEV_TO_PI * args.emin)
pi_max = int(KEV_TO_PI * args.emax)
cmd  = 'fselect %s %s ' % (args.inputevtfits,args.outputevtfits)
cmd += '"(PI >= %d).and.(PI < %d)" ' % (pi_min,pi_max)
print(cmd);os.system(cmd)
