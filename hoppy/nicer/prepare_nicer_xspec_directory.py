#!/usr/bin/env python

import os 
import sys 
import glob
import argparse

#import hoppy.nicer.nipipeline as nip

if __name__=="__main__":

	sys.stdout.write('... prepare xspec fitting directories ...\n')

	parser = argparse.ArgumentParser(
		prog='prepare_nicer_xspec_directory.py',
		usage='prepare_nicer_xspec_directory.py obsid_path_str [--outdir (default=nispec)] [--flag_evt]',
		description='prepare xspec fitting directory',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'obsid_path_str',metavar='obsid_path_str',type=str,        
		help='Input ObsID directory for glob. e.g., "../../../pipe/v180809_ws1/niout/o1u20cut/10500501*"')
	parser.add_argument(
		'--outdir',metavar='--outdir',type=str,default="nispec",
		help='Output directory (default=nispec).')
	parser.add_argument('--flag_evt', action='store_true', 
		help='flag to copy of event files.') 	
	args = parser.parse_args()	

	for obsid_path in glob.glob(args.obsid_path_str):
		print(obsid_path)
		obsid = os.path.basename(obsid_path)
		outdir = '%s/%s' % (args.outdir,obsid)
		cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
		print(cmd);os.system(cmd)

		indir = '%s/work/*/process/speclc/product/' % obsid_path
		for pha in glob.glob('%s/*.pha' % indir):
			cmd = 'cp %s %s' % (pha,outdir)
			print(cmd);os.system(cmd)
		for pdf in glob.glob('%s/*.pdf' % indir):
			cmd = 'cp %s %s' % (pdf,outdir)
			print(cmd);os.system(cmd)
		if args.flag_evt:
			print("evt")
			for evt in glob.glob('%s/work/o1u20cut/process/*evt.gz' % obsid_path):			
				cmd = 'cp %s %s' % (evt,outdir)
				print(cmd);os.system(cmd)
