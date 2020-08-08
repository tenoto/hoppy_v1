#!/usr/bin/env python

import os
import argparse
import numpy as np 

from astropy.io import fits
from astropy.time import Time

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-08 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('nibgdlc3c50.py',
#		usage='%(prog)s ',
		description="""
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--tbin', type=float, default=None,
		help='tbin')		
	parser.add_argument('--bkgidxdir', type=str, default=None,
		help='bkgidxdir (see fhelp nibackgen3C50)')	
	parser.add_argument('--bkglibdir', type=str, default=None,
		help='bkglibdir (see fhelp nibackgen3C50)')	
	parser.add_argument('--gainepoch', type=str, default=None,
		help='gainepoch (see fhelp nibackgen3C50)')		
	parser.add_argument('--ufafile', type=str, default=None,
		help='ufafile (see fhelp nibackgen3C50)')	
	parser.add_argument('--clfile', type=str, default=None,
		help='clfile (see fhelp nibackgen3C50)')			
	parser.add_argument('--dtmin', type=float, default=None,
		help='dtmin (see fhelp nibackgen3C50)')	
	parser.add_argument('--dtmax', type=float, default=None,
		help='dtmax (see fhelp nibackgen3C50)')	
	parser.add_argument('--hbgcut', type=float, default=None,
		help='hbgcut (see fhelp nibackgen3C50)')	
	parser.add_argument('--fpmofflist', type=str, default=None,
		help='fpmofflist (see fhelp nibackgen3C50)')											
	return parser

def run(args):
	print(args)

	clevt = args.clfile
	clhdu = fits.open(clevt)

	tstart = clhdu['EVENTS'].header['TSTART']
	tstop = clhdu['EVENTS'].header['TSTOP']

	tcenter_list = []
	i = 0
	tcenter = tstart 
	while tcenter < tstop:
		tcenter = tstart + 0.5 * args.tbin + args.tbin * i
		tcenter_list.append(tcenter)
		i += 1
	tcenter = np.array(tcenter_list)

	tcenter_gti_list = []
	for t in tcenter:
		for gti in clhdu['GTI'].data:
			if t > gti[0] and t <= gti[1]:
				tcenter_gti_list.append(t)
	gti_tcenter = np.array(tcenter_gti_list)

	outdir = 'nibgdlc3c50'
	for i in range(len(gti_tcenter)):
		print(i)
		gti_start = gti_tcenter[i] - 0.5 * args.tbin
		gti_stop = gti_tcenter[i] + 0.5 * args.tbin		
		gtinum = i + 1
		subdir = '%s/gti%03d' % (outdir,gtinum)
		cmd = 'rm -rf %s;\n' % subdir
		cmd += 'mkdir -p %s' % subdir
		print(cmd);os.system(cmd)

		fgti_txt = '%s/gti%03d.txt' % (subdir,gtinum)
		f = open(fgti_txt,'w')
		f.write('%.10f	%.10f\n' % (gti_start,gti_stop))
		f.close()

		fgti_fits = '%s/gti%03d.fits' % (subdir,gtinum)
		cmd = 'fconv_txt2gti.py -i %s -o %s ' % (fgti_txt,fgti_fits)
		print(cmd);os.system(cmd)
		
def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	run(args)

if __name__=="__main__":
	main()