#!/usr/bin/env python

import os 
import sys 
import glob
import astropy.io.fits as pyfits 
import argparse

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='run_photonphase_multi_observations.py',
		usage='run_photonphase_multi_observations.py obsid_path_str parfile [--outdir (default=nipulse)]',
		description='run series of photonphase for NICER data.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'obsid_path_str',metavar='obsid_path_str',type=str,        
		help='Input ObsID directory for glob. e.g., "../../../pipe/v180809_ws1/niout/o1u20cut/10500501*"')
	parser.add_argument(
		'parfile',metavar='parfile',type=str,        
		help='Input folding parameter file (TEMPO format).')	
	parser.add_argument(
		'--filtername',metavar='--filtername',type=str,default="o1u20cut",
		help='filter verion directory name, e.g., data/1020130105/work/[o1u20cut] (default=o1u20cut).')	
	parser.add_argument(
		'--outdir',metavar='--outdir',type=str,default="nipulse",
		help='Output directory (default=nipulse).')
	args = parser.parse_args()	

	print(args)

	for obsid_path in glob.glob(args.obsid_path_str):
		obsid = os.path.basename(obsid_path)
		print("processing %s" % obsid)
		outdir = '%s/%s' % (args.outdir,obsid)
		
		inevt_list = glob.glob('data/%s/work/%s/process/ni%s_0mpu7_cl_%s.evt.gz' % (obsid,args.filtername,obsid,args.filtername))
		if len(inevt_list) == 0:
			sys.stderr.write('file does not exists.\n')
			continue
		inevt = inevt_list[0]
		hdu = pyfits.open(inevt)
		if len(hdu[1].data) == 0:
			sys.stderr.write('no events.\n')
			continue
		cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
		print(cmd);os.system(cmd)

		orbfile = '%s/auxil/ni%s.orb.gz' % (obsid_path,obsid)
		outevt  = '%s/ni%s_0mpu7_cl_%s_photonphase.evt' % (outdir,obsid,args.filtername)
		plotfile = outevt.replace('.evt','.pdf')
		logfile =  outevt.replace('.evt','.log')

		fname_script = outevt.replace('.evt','.sh')
		f = open(fname_script,'w')
		f.write('#!/bin/sh -f\n')
		dump  = 'photonphase '
		dump += '--orbfile %s ' % orbfile
		dump += '--plotfile %s ' % plotfile
		dump += '--addphase '
		dump += '--outfile %s ' % outevt
		dump += '%s ' % inevt
		dump += '%s  ' % args.parfile
		dump += '> %s 2>&1' % logfile 
		f.write(dump)
		f.close()
		cmd = 'chmod +x %s' % fname_script
		print(cmd);os.system(cmd)

		print("...now photonphase is running for %s" % obsid)
		cmd = './%s' % fname_script
		print(cmd);os.system(cmd)		






