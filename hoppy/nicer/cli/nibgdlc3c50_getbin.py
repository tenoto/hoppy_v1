#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np 

from astropy.io import fits
#from astropy.time import Time

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-17 : original version

def get_parser():
	parser = argparse.ArgumentParser('nibgdlc3c50_getbin.py',
		description="""
Generate a csv format file including a time series of tstart and tstop for nibgdlc3c50.py 
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--clfile', type=str, default=None, required=True,
		help='cleaned evetn file for GTIs')		
	parser.add_argument('--tbin', type=float, default=None, required=True,
		help='time bin (sec)')				
	parser.add_argument('--outfile', type=str, default=None, required=True,
		help='output csv format file for tstart and tstop series')		
	return parser

def nibgdlc3c50_getbin(args):
	print("\n %s" % (sys._getframe().f_code.co_name))

	print(args)

	outdir = 'tmp_nibgdlc3c50_getbin'
	cmd = 'rm -rf %s;mkdir -p %s' % (outdir,outdir)
	print(cmd);os.system(cmd)

	clevt = args.clfile
	clhdu = fits.open(clevt)

	tbin = float(args.tbin)
	tstart = clhdu['EVENTS'].header['TSTART']
	tstop  = clhdu['EVENTS'].header['TSTOP']

	f = open(args.outfile,'w')
	dump = 'START_MET,STOP_MET,Exposure_sec\n'
	f.write(dump)

	gti_tstart = tstart
	gti_tstop = tstart
	i = 0
	while gti_tstop < tstop + tbin:
		gti_tstart = tstart + tbin * i 
		gti_tstop = tstart + tbin * (i+1)

		flag_gti_in = False
		for gti in clhdu['GTI'].data:
			if gti['START'] < gti_tstart < gti['STOP'] or gti['START'] < gti_tstop < gti['STOP']:
				flag_gti_in = True
				break 
		if not flag_gti_in:
			i += 1 
			continue

		fgti_txt = '%s/tmp_gti%05d.txt' % (outdir,i)
		ftmp = open(fgti_txt,'w')
		ftmp.write('%.10f	%.10f\n' % (gti_tstart,gti_tstop))
		ftmp.close()

		fgti_fits = '%s/tmp_gti%05d.fits' % (outdir,i)
		cmd = 'fconv_txt2gti.py -i %s -o %s ' % (fgti_txt,fgti_fits)
		print(cmd);os.system(cmd)
		
		fgti_merge = '%s/tmp_gti%05d_merge.fits' % (outdir,i)
		cmd = 'mgtime "%s,%s+2" %s AND' % (fgti_fits,args.clfile,fgti_merge)
		print(cmd);os.system(cmd)

		mghdu = fits.open(fgti_merge)
		gti_exposure = sum(mghdu['STDGTI'].data['STOP']-mghdu['STDGTI'].data['START'])
		f.write('%d,%d,%.1f\n' % (gti_tstart, gti_tstop,gti_exposure))
		i += 1	

	f.close()

	cmd = 'rm -rf %s' % (outdir)
	print(cmd);os.system(cmd)	

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	nibgdlc3c50_getbin(args)

if __name__=="__main__":
	main()