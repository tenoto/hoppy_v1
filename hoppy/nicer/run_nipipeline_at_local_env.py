#!/usr/bin/env python

import os 
import sys 
import glob
import argparse

import hoppy.nicer.nipipeline as nip

if __name__=="__main__":

	sys.stdout.write('... run a single [multiple] ObsID process[es] at local environemnt ...\n')

	parser = argparse.ArgumentParser(
		prog='run_nipipeline_at_local_env.py',
		usage='run_nipipeline_at_local_env.py obsid [fyaml_param.yaml] [outdir] [flag_recreate]',
		description='A series of interfaces for NICER pipelines.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'obsid',metavar='obsid',type=str,        
		help='Input ObsID. e.g., "1012010133" for a single or "10120101*" for multiple data.')
	parser.add_argument(
		'--fyaml_param',metavar='fyaml_param',type=str,     
		default=os.getenv('NICER_DEFAULT_FPARAM_FILE'),   
		help='yaml file for input parameters.')
	parser.add_argument(
		'-r', '--recreate',action='store_true',dest='flag_recreate',
		default=False,
		help='recreate flag.')
	args = parser.parse_args()	

	obsid_path_list = glob.glob('%s/*/%s' % (os.getenv('NICER_DATA_REPOSITORY'),args.obsid))		
	print(obsid_path_list)
	
	for obsid_path in obsid_path_list:
		#cmd  = 'nipipeline.py %s ' % obsid_path
		#cmd += '--fyaml_param %s ' % args.fyaml_param
		#if args.recreate:
		#	cmd += '--recreate '
		#print(cmd)
		niobsid = nip.NiObsID(obsid_path,args.fyaml_param)
		niobsid.run(args.flag_recreate)