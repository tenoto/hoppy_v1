#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run following commands (similar to nicerl2)
1. niprefilter2 - create augmented NICER-specific filter file
2. nicercal - apply standard NICER calibration
3. nimaketime - create standard screening good time intervals
4. nicermergeclean - combine per-MPU data and filter/screen
"""

__author__  = 'Teru Enoto'
__date__    = '2018 August 9'
__version__ = '1.00'
__keyword__ = 'nipipeline_v%s' % __version__

import os
import sys  
import yaml 
import glob 
import argparse

NICER_SOFT_PATH = '%s/hoppy/nicer' % os.getenv('HOPPY_PATH')
DEFAULT_FPARAM_FILE = '%s/nipipeline/pipeline_setup_default.yaml' % NICER_SOFT_PATH

class NiObsID():
	def __init__(self,obsid_path,fyaml_param):
		self.obsid_path = obsid_path
		self.fyaml_param = fyaml_param
		self.obsid = os.path.basename(self.obsid_path)

	def check_input_files(self):
		if not os.path.exists(self.obsid_path):
			sys.stderr.write('[error] obsid_path does not exist: %s\n' % self.obsid_path)			
			quit()				
		if not os.path.exists(self.fyaml_param)	:
			sys.stderr.write('[error] fyaml_param does not exist: %s\n' % self.fyaml_param)			
			quit()				
		try:
			self.param = yaml.load(open(self.fyaml_param))
		except:
			sys.stderr.write('[error] fyaml_param can not be loaded.: %s\n' % self.fyaml_param)			
			quit()				
		self.outdir = self.param['outdir']
		self.outdir_obsid = '%s/%s' % (self.outdir,self.obsid)

		sys.stdout.write('===============================\n')
		sys.stdout.write('NiObsID\n')		
		sys.stdout.write('obsid_path: %s\n' % self.obsid_path)				
		sys.stdout.write('outdir: %s\n' % self.outdir)										
		sys.stdout.write('fyaml_param: %s\n' % self.fyaml_param)						
		print(self.param)		
		sys.stdout.write('===============================\n')		

#	def prepare_directory(self):
#		obsiddir = '%s/%s' % (self.outdir,obsid)

	def run(self):
		self.check_input_files()


if args.flag_recreate:
	cmd = 'rm -rf %s' % (obsiddir)
	print(cmd);os.system(cmd)

if not os.path.exists(obsiddir):
	cmd = 'mkdir -p %s/xti' % (obsiddir)
	print(cmd);os.system(cmd)

	pwd = os.getcwd()
	os.chdir(obsiddir)
	cmd  = 'cp -r %s/auxil .;' % indir_path
	cmd += 'ln -s %s/log .;' % indir_path
	print(cmd);os.system(cmd)
	os.chdir('xti')
	cmd  = 'ln -s %s/xti/hk .;' % indir_path
	cmd += 'ln -s %s/xti/event_uf .;' % indir_path
	cmd += 'mkdir event_cl;' 
	#cmd += 'cp -r %s/xti/event_cl .;' % indir_path
	print(cmd);os.system(cmd)
	os.chdir(pwd)

	def run(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

if __name__=="__main__":

	sys.stdout.write('... run a single ObsID process (%s) ...\n' % __keyword__)
	sys.stdout.write('DEFAULT_FPARAM_FILE: %s\n' % DEFAULT_FPARAM_FILE)

	parser = argparse.ArgumentParser(
		prog='nipipeline.py',
		usage='nipipeline.py obsid_path [fyaml_param.yaml] [outdir] [flag_recreate]',
		description='A series of interfaces for NICER pipelines.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'obsid_path',metavar='obsid_path',type=str,        
		help='Input ObsID directory path. The directory should be a single NICER observation directory, which in turn contains xti/{events_uf,events_cl,hk,auxil} subdirectories.')
	parser.add_argument(
		'--fyaml_param',metavar='fyaml_param',type=str,     
		default=DEFAULT_FPARAM_FILE,   
		help='yaml file for input parameters.')
	parser.add_argument(
		'--outdir', metavar='outdir',type=str,     
		default="niout",
		help='main output directory.')
	parser.add_argument(
		'-r', '--recreate',action='store_true',dest='flag_recreate',
		default=False,
		help='recreate flag.')	
	args = parser.parse_args()	
	
	niobsid = NiObsID(args.obsid_path,args.fyaml_param)
	niobsid.run()



