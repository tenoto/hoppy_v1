#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# HISTORY
# 2018-09-19: Modified to new version (Teru Enoto)
# 2018-05-08: Generated from scratch (Teru Enoto)

__author__  = 'Teru Enoto'
__date__    = '2018 September 19'
__version__ = '1.01'

import os 
import sys 
import argparse
import numpy as np 
import astropy.io.fits as pyfits

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rc
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['font.size'] = 18
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.xmargin'] = '0' #'.05'
plt.rcParams['axes.ymargin'] = '.2'

def generate_cmap(colors):
	"""
	https://qiita.com/kenmatsu4/items/fe8a2f1c34c8d5676df8
	"""
	values = range(len(colors))
	vmax = np.ceil(np.max(values))
	color_list = []
	for v, c in zip(values, colors):
		color_list.append( ( v/ vmax, c) )
	return LinearSegmentedColormap.from_list('custom_cmap', color_list)

def fplot_hist2d_pulseprofile(infits,outfits,
	emin=0.1,emax=10.0,ebin=20,pbin=20,pcol="PULSE_PHASE"):
	sys.stdout.write("---fplot_hist2d_pulseprofile---\n")
	sys.stdout.write("infits: %s \n" % infits)
	sys.stdout.write("outfits: %s \n" % outfits)
	sys.stdout.write("emin: %.3f (keV)\n" % emin)
	sys.stdout.write("emax: %.3f (keV)\n" % emax)	
	sys.stdout.write("ebin: %d\n" % ebin)
	sys.stdout.write("pbin: %d\n" % pbin)
	sys.stdout.write("pcol: %s\n" % pcol)

	if not os.path.exists(infits):
		sys.stdout.write('Error: file does not exist. %s' % infits)
		quit()		

	outpdf = outfits.replace('.fits','.pdf')

	hdu = pyfits.open(infits)
	nevents = len(hdu['EVENTS'].data)
	if nevents == 0:
		sys.stdout.write('Skip: No events in fits file %s' % infits)
		quit()
	telescope  = hdu[0].header['TELESCOP']
	instrument = hdu[0].header['INSTRUME']
	KEV_TO_PI = None
	if telescope == 'NICER' and instrument == 'XTI':
		KEV_TO_PI = 1000.0/10.0 # Conversion from PI to keV (PI is in units of 10 eV)
	PI_TO_KEV = 1/KEV_TO_PI
	sys.stdout.write('TELESCOP: %s\n' % telescope)
	sys.stdout.write('INSTRUME: %s\n' % instrument)
	sys.stdout.write('KEV_TO_PI: %.3f\n' % KEV_TO_PI) 
	sys.stdout.write('PI_TO_KEV: %.3f\n' % PI_TO_KEV) 	
	sys.stdout.write('NEVENTS: %d\n' % nevents)

	phase  = hdu[1].data['PULSE_PHASE']
	energy = (hdu[1].data['PI'] + np.random.rand(nevents)) * PI_TO_KEV
	x = np.append(phase, phase+1.0)
	y = np.append(energy, energy)

	# https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist2d.html
	fig = plt.figure(figsize=(10,8))
	#cm = generate_cmap(['#00008B', '#aaaaab', '#FFFFFF', '#F4D793', '#F4A460'])
	hist2d, xedges, yedges, image = plt.hist2d(x,y,
		#norm=mpl.colors.PowerNorm(gamma=1./2.),
		norm=mpl.colors.LogNorm(),
		cmin=1.0,#vmin=4000,
		#cmap = cm,
		cmap='inferno',
		bins=[pbin,ebin],range=[[0.0,2.0],[emin,emax]])
	plt.xlabel('Pulse Phase')
	plt.ylabel('Energy (keV)')
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Photon Counts')
	fig.savefig(outpdf)

if __name__=="__main__":

	usage = """
Two-dimensional plotting 	
"""

	parser = argparse.ArgumentParser(
		prog='fplot_hist2d_pulseprofile.py',
		usage='fplot_hist2d_pulseprofile.py infits outfits [--emin] [--emax] [--ebin] [--pbin] [--pcol]',
		description=usage,
		epilog='',
		add_help=True)

	parser.add_argument('infits',metavar='infits',type=str,        
		help='Input event fits file. Starting @ shows a list.')	
	parser.add_argument('outfits',metavar='outfits',type=str, 		       
		help='Output 2d histogram fits file.')	
	parser.add_argument('--emin','--energy_min',metavar='emin',
		type=float,default=0.0,
		help='Energy min (keV) [default=0.1 keV].')		
	parser.add_argument('--emax','--energy_max',metavar='emax',
		type=float,default=10.0,
		help='Energy max (keV) [default=10.0 keV].')
	parser.add_argument('--ebin','--energy_nbin',metavar='ebin',
		type=int,default=20,
		help='Nbin for the energy direction [default=20].')		
	parser.add_argument('--pbin','--phase_nbin',metavar='pbin',
		type=int,default=20,
		help='Nbin for the phase direction [default=20].')	
	parser.add_argument('--pcol','--phase_column',metavar='pcol',
		type=str,default="PULSE_PHASE",
		help='Phase column name [default=PULSE_PHASE].')											
	args = parser.parse_args()	
	print(args)

	fplot_hist2d_pulseprofile(
		infits=args.infits,
		outfits=args.outfits,
		emin=args.emin,emax=args.emax,ebin=args.ebin,
		pbin=args.pbin,pcol=args.pcol)




	
