#!/usr/bin/env python

__author__  = 'Teruaki Enoto'
__version__ = '0.01'
__date__    = '2018 September 28'

import os 
import sys 
import astropy.io.fits as fits
import argparse

import pandas as pd 
import matplotlib as mpl
import matplotlib.pyplot as plt 

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '14'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.05' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

def compare_gti(fname_input_file_lst,outfname,
	start_colname='BARY_START',stop_colname='BARY_STOP',xlabel='Barycentric MJD'):
	input_file_list = []
	observatory_list = []
	for line in open(fname_input_file_lst):
		infile = line.split(',')[0]
		if not os.path.exists(infile):
			sys.stderr.write('infile %s does not exist.' % infile)
			quit()
		input_file_list.append(infile)
		observatory_list.append(line.split(',')[1].strip())

	fig, axes = plt.subplots(1,1,figsize=(9.6,3.6))

	mjd_min = 1e+10 
	mjd_max = 0.0
	for i in range(len(input_file_list)):
		hdu = fits.open(input_file_list[i])
		for line in hdu['GTI'].data:
			if line[start_colname] < mjd_min:
				mjd_min = line[start_colname]
			if line[stop_colname] > mjd_max:
				mjd_max = line[stop_colname]
	mjd_width = mjd_max - mjd_min
	mjd_offset = int(mjd_min)

	for i in range(len(input_file_list)):
		hdu = fits.open(input_file_list[i])
		for line in hdu['GTI'].data:
			plt.hlines(i,line[start_colname]-mjd_offset,line[stop_colname]-mjd_offset,
				colors='C%d' % i,label=observatory_list[i],lw=5)
		plt.text(mjd_min-mjd_offset-0.18*mjd_width,i-0.1,observatory_list[i])

	plt.ylim(-1.0,len(input_file_list))
	plt.xlabel('%s - %.1f' % (xlabel,mjd_offset) )
	plt.yticks([])
	plt.savefig(outfname)

if __name__=="__main__":

	usage = """
This script compare GTIs and plot them.
"""

	parser = argparse.ArgumentParser(
		prog='fcompare_gti.py',
		usage='fcompare_gti.py input_file.lst',
		description=usage,epilog='',add_help=True)
	parser.add_argument('fname_input_file_lst',metavar='fname_input_file_lst',type=str,help='Input GTI fits files.')
	parser.add_argument('outfname',metavar='outfname',type=str,help='Output file name.')	
	parser.add_argument('--start_colname',metavar='start_colname',
		type=str,help='START column name.',default='BARY_START')	
	parser.add_argument('--stop_colname',metavar='stop_colname',
		type=str,help='START column name.',default='BARY_STOP')		
	args = parser.parse_args()	

	compare_gti(args.fname_input_file_lst,args.outfname,
		start_colname=args.start_colname,
		stop_colname=args.stop_colname
		)

