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
__version__ = '1.01'
__keyword__ = 'nipipeline_v%s' % __version__

import os
import sys  
import yaml 
import glob 
import argparse
import astropy.io.fits as pyfits

NICER_SOFT_PATH = '%s/hoppy/nicer' % os.getenv('HOPPY_PATH')
DEFAULT_FPARAM_FILE = '%s/nipipeline/pipeline_setup_default.yaml' % NICER_SOFT_PATH

class NiObsID():
	def __init__(self,obsid_path,fyaml_param):
		self.obsid_path = obsid_path
		self.fyaml_param = fyaml_param
		self.obsid = os.path.basename(self.obsid_path)

	def check_input_files(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)		

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

	def prepare_directory(self,flag_recreate=False):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)		

		if flag_recreate:
			cmd = 'rm -rf %s' % (self.outdir_obsid)
			print(cmd);os.system(cmd)

		if not os.path.exists(self.outdir_obsid):
			cmd = 'mkdir -p %s/xti' % (self.outdir_obsid)
			print(cmd);os.system(cmd)

		pwd = os.getcwd()
		os.chdir(self.outdir_obsid)
		cmd  = 'cp -r %s/auxil .;' % self.obsid_path
		cmd += 'ln -s %s/log .;' % self.obsid_path
		print(cmd);os.system(cmd)
		os.chdir('xti')
		cmd  = 'ln -s %s/xti/hk .;' % self.obsid_path
		cmd += 'ln -s %s/xti/event_uf .;' % self.obsid_path
		cmd += 'mkdir event_cl;' 
		#cmd += 'cp -r %s/xti/event_cl .;' % indir_path
		print(cmd);os.system(cmd)
		os.chdir(pwd)

	def run_nicercal(self):
		"""
		1. nicercal - apply standard NICER calibration
		"""
		self.fname_mkffile = '%s/auxil/ni%s.mkf.gz' % (self.outdir_obsid,self.obsid)

		self.fname_log_nicercal = '%s/xti/event_cl/1_nicercal.log' % self.outdir_obsid
		self.fname_outfilefile = '%s/xti/event_cl/1_nicercal_ufalist.lis' % self.outdir_obsid
		cmd = 'nicercal indir=%s outfilefile=%s >& %s' % (self.outdir_obsid,self.fname_outfilefile,self.fname_log_nicercal)
		print(cmd);
		if os.path.exists(self.fname_log_nicercal):
			print("...already processed. skipped.")
		else:
			os.system(cmd)
			cmd = 'gzip %s' % self.fname_mkffile.replace('.gz','')
			print(cmd);os.system(cmd)

		self.fname_ufamerge_evt = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt' % (self.outdir_obsid,self.obsid)
		cmd = 'nimpumerge infiles=@%s outfile=%s ' % (self.fname_outfilefile,self.fname_ufamerge_evt)
		print(cmd);
		if os.path.exists(self.fname_ufamerge_evt+'.gz'):
			print("...already processed. skipped.")
		else:
			os.system(cmd)
			cmd = 'gzip %s' % self.fname_ufamerge_evt
			print(cmd);os.system(cmd)
		self.fname_ufamerge_evt += '.gz'

	def run_niprefilter2(self):
		"""
		2. niprefilter2 - create augmented NICER-specific filter file
		"""
		self.fname_log_niprefilter2 = '%s/auxil/2_niprefilter2.log' % self.outdir_obsid

		if os.path.exists(self.fname_log_niprefilter2):
			print("...niprefilter2. already processed. skipped.")
		else:
			cmd = 'gunzip %s' % self.fname_mkffile
			print(cmd);os.system(cmd)

			cmd = 'niprefilter2 clobber=yes indir=%s infile=%s outfile=INFILE >& %s' % (self.outdir_obsid,self.fname_mkffile.replace('.gz',''),self.fname_log_niprefilter2)
			print(cmd);os.system(cmd)

			cmd = 'gzip %s' % self.fname_mkffile.replace('.gz','')
			print(cmd);os.system(cmd)

	def run_nimaketime(self):
		"""
		3. nimaketime - create standard screening good time intervals
		"""
		self.dir_work  = '%s/work/%s' % (self.outdir_obsid,self.param['gtiexpr_basestr'])
		self.dir_proc = '%s/process' % self.dir_work
		self.dir_log  = '%s/log' % self.dir_proc
		cmd = 'mkdir -p %s' % self.dir_log
		print(cmd);os.system(cmd)

		cmd = 'cp %s %s' % (self.fyaml_param,self.dir_log)
		print(cmd);os.system(cmd)

		self.fname_log_nimaketime = '%s/3_nimaketime.log' % self.dir_log
		self.fname_expr_nimaketime = '%s/3_nimaketime.expr' % self.dir_log
		self.fname_gtifile = '%s/ni%s_%s.gti' % (self.dir_proc,self.obsid,self.param['gtiexpr_basestr'])
		cmd = 'nimaketime infile=%s outfile=%s expr="%s" outexprfile="%s" >& %s' % (self.fname_mkffile, 
			self.fname_gtifile, self.param['nimaketime_expr'], 
			self.fname_expr_nimaketime, self.fname_log_nimaketime)
		print(cmd);os.system(cmd)

	def run_nicermergeclean(self):
		"""
		4. nicermergeclean - combine per-MPU data and filter/screen
		"""
		self.fname_log_nicermergeclean = '%s/4_nicermergeclean.log' % self.dir_log
		self.fname_ufaevt = '%s/ni%s_0mpu7_ufa_%s.evt' % (self.dir_proc,self.obsid,self.param['gtiexpr_basestr'])
		self.fname_clevt  = '%s/ni%s_0mpu7_cl_%s.evt' % (self.dir_proc,self.obsid,self.param['gtiexpr_basestr'])
		#self.fname_ufaevt = '%s/xti/event_cl/ni%s_0mpu7_ufa_%s.evt' % (dir_work,obsid,self.param['gtiexpr_basestr'])
		#self.fname_clevt  = '%s/xti/event_cl/ni%s_0mpu7_cl_%s.evt' % (dir_work,obsid,self.param['gtiexpr_basestr'])
		cmd  = 'nicermergeclean infiles=@%s ' % self.fname_outfilefile
		cmd += 'ufafile=%s ' % self.fname_ufaevt
		cmd += 'clfile=%s ' % self.fname_clevt
		cmd += 'gtifile=%s ' % self.fname_gtifile
		cmd += '>& %s ' % self.fname_log_nicermergeclean
		print(cmd);os.system(cmd)

		cmd  = 'fparkey %s %s NIPIPVER comm="nipipeline.py version" add=yes;'  % (__version__,self.fname_ufaevt)
		cmd += 'fparkey %s %s NIPIPVER comm="nipipeline.py version" add=yes;'  % (__version__,self.fname_clevt)
		print(cmd);os.system(cmd)

		cmd = 'gzip %s;\n' % self.fname_ufaevt
		cmd += 'gzip %s;\n' % self.fname_clevt
		print(cmd);os.system(cmd)

		self.fname_ufaevt += '.gz'
		self.fname_clevt += '.gz'

	def run_make_bgdspec(self):
		"""
		5. make_bgdspec_mit3C50.py
		"""
		self.dir_speclc = '%s/speclc' % self.dir_proc
		prefix = 'ni%s_%s' % (self.obsid,self.param['gtiexpr_basestr'])

		self.fname_log_nibgdspec_mit3C50 = '%s/5_nibgdspec_mit3C50.log' % self.dir_log

		#cmd  = 'make_bgdspec_mit3C50.py '
		cmd  = 'nibgdspec_mit3C50.py '
		cmd += '%s %s ' % (self.fname_ufaevt,self.fname_clevt)
		cmd += '--outdir %s ' % self.dir_speclc
		cmd += '--prefix %s ' % prefix
		cmd += '--exclude %s ' % str(self.param['exclude_detid']).replace('[','').replace(']','').replace(' ','')
		cmd += '--tbin %.1f ' % self.param['mitbgd_tbin']
		cmd += '--lctbin %.1f --lcemin %.1f --lcemax %.1f ' % (
        	self.param['mitbgd_ql_lc_tbin'],self.param['mitbgd_ql_lc_emin'],self.param['mitbgd_ql_lc_emax'])
		cmd += '> %s ' % (self.fname_log_nibgdspec_mit3C50)
		print(cmd);os.system(cmd)

	def plot_energy_selected_curve(self):
		"""
		6. Plot energy-selected curve
		"""
		if len(glob.glob('%s/product' % self.dir_speclc)) == 0:
			sys.stderr.write("no gti. exit.\n")
			exit()
		self.product_cl_evt = glob.glob('%s/product/ni%s_*_clscr_gtisel.evt' % (self.dir_speclc,self.obsid))[0]

		self.flc_eband_list = []
		for eband in self.param['lc_eband_list']:
			fname_flc_eband = '%s_%sto%skeV.flc' % (self.product_cl_evt.replace('.evt',''),str(eband[0]).replace('.','p'),str(eband[1]).replace('.','p'))
			cmd  = 'fselect_filter_energy.py '
			cmd += '%s ' % self.product_cl_evt # input event fits file
			cmd += '%s ' % fname_flc_eband # output event fits file
			cmd += '%.2f ' % eband[0] # energy min (keV)
			cmd += '%.2f ' % eband[1] # energy max (keV)
			print(cmd);os.system(cmd)
			self.flc_eband_list.append(fname_flc_eband)

		self.product_bandflc = self.product_cl_evt.replace('.evt','_eband.flc')
		self.product_bandflc_ps = self.product_bandflc.replace('.flc','_flc.ps')
		cmd = 'lcurve %d <<EOF\n' % len(self.flc_eband_list)
		for flc_name in self.flc_eband_list:
			cmd += '%s\n' % flc_name
		cmd += '-\n'	
		cmd += '%d\n' % self.param['mitbgd_ql_lc_tbin']
		cmd += '%d\n' % self.param['lc_nbint']
		cmd += '%s\n' % self.product_bandflc
		cmd += 'yes\n'
		cmd += '/xw\n'
		cmd += '%d\n' % len(self.flc_eband_list)
		i = 2
		for eband in self.param['lc_eband_list']:
			cmd += 'lab y%d %.1f-%.1f keV\n' % (i,eband[0],eband[1])
			cmd += 'col %d on %d\n' % (i,i)
			cmd += 'lab rotate\n'
			i += 1 
		cmd += 'la f %s (%s)\n' % (self.obsid,self.param['gtiexpr_basestr'])
		cmd += 'lwid 5 on 1..100\n'
		cmd += 'lab rotate\n'
		cmd += 'hard %s/cps\n'	% self.product_bandflc_ps
		cmd += 'exit'
		print(cmd);os.system(cmd)
		cmd = 'ps2pdf %s' % self.product_bandflc_ps
		print(cmd);os.system(cmd)

		cmd = 'rm -f %s; mv %s %s' % (self.product_bandflc_ps,
			os.path.basename(self.product_bandflc_ps).replace('.ps','.pdf'),
			os.path.dirname(self.product_bandflc_ps))
		print(cmd);os.system(cmd)		

	def calculate_count_rates(self):
		"""
		7. Calculate Counts 
		"""
		hdu = pyfits.open(self.product_cl_evt)
		target  = hdu[0].header['OBJECT']
		dateobs = hdu[0].header['DATE-OBS']
		exposure = hdu[1].header['EXPOSURE']
		title = '%s %s %s (%.1f s)' % (target,self.obsid,dateobs,exposure)

		self.fname_read_xcm = glob.glob('%s/product/ni%s_*_read.xcm' % (self.dir_speclc,self.obsid))[0]
		self.fname_residual_ps = self.product_cl_evt.replace('.evt','_spec_residual.ps')
		cmd  = 'xspec<<EOF\n'
		cmd += 'setplot device tmp/null\n'			
		cmd += '@%s\n' % self.fname_read_xcm
		cmd += 'setplot energy\n'
		cmd += 'notice **-**\n'
		cmd += 'ignore **-0.2 15.0-**\n'
		cmd += 'setplot rebin 3 50 1 \n'
		cmd += 'setplot rebin 3 50 2 \n'
		cmd += 'ipl d\n'
		cmd += 'r y -3.0 3.0\n'
		cmd += 'r x 0.2 15.0\n'
		cmd += 'lab 1 pos 0.2 0 loc 0 1 " " ls 3 col 4 \n'
		cmd += 'lwid 5 \n'
		cmd += 'lwid 5 on 1..100\n'
		cmd += 'time off\n'
		cmd += 'la t %s\n' % title 
		cmd += 'la f %s (red=obs, green=bgd_model, black=bgd_subtracted)\n' % self.obsid 
		cmd += 'hard %s/cps\n' % self.fname_residual_ps
		cmd += 'exit\n'
		cmd += 'data 2 none\n'
		tmp_logfile_list = []
		for eband in self.param['xspec_eband_list']:
			logfile = 'tmp_rate_%sto%skeV.log' % (str(eband[0]).replace('.','p'),str(eband[1]).replace('.','p'))
			cmd += 'notice **-**\n'
			cmd += 'log %s\n' % logfile
			cmd += 'ignore **-%.2f %.2f-**\n' % (eband[0],eband[1])
			cmd += 'show rate\n'
			cmd += 'log none\n'
			tmp_logfile_list.append(logfile)
		cmd += 'exit\n'	
		print(cmd);os.system(cmd)

		cmd = 'ps2pdf %s' % self.fname_residual_ps
		print(cmd);os.system(cmd)

		cmd = 'rm -f %s; mv %s %s' % (self.fname_residual_ps,
			os.path.basename(self.fname_residual_ps).replace('.ps','.pdf'),
			os.path.dirname(self.fname_residual_ps))
		print(cmd);os.system(cmd)

		self.xspec_rate_werror_list = []
		for logfile in tmp_logfile_list:
			i = tmp_logfile_list.index(logfile)
			emin = self.param['xspec_eband_list'][i][0]
			emax = self.param['xspec_eband_list'][i][1]
			for line in open(logfile):
				cols = line.split()
				if cols[0] == '#Net':
					rate = float(cols[6])
					rate_error = float(cols[8])
					print(emin,emax,rate,rate_error)
					self.xspec_rate_werror_list.append([rate,rate_error])
			cmd = 'rm -f %s\n' % logfile
			print(cmd);os.system(cmd)

		self.param['xspec_rate_werror_list'] = self.xspec_rate_werror_list
		self.fname_result_yaml = self.fname_read_xcm.replace('_read.xcm','_out.yaml')
		with open(self.fname_result_yaml, 'w') as yaml_file:
		    yaml.dump(self.param, yaml_file, default_flow_style=True)

	def run(self,flag_recreate):
		self.check_input_files()
		self.prepare_directory(flag_recreate)
		self.run_nicercal()
		self.run_niprefilter2()
		self.run_nimaketime()
		self.run_nicermergeclean()
		self.run_make_bgdspec()
		self.plot_energy_selected_curve()
		self.calculate_count_rates()

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
		'-r', '--recreate',action='store_true',dest='flag_recreate',
		default=False,
		help='recreate flag.')	
	args = parser.parse_args()	
	
	niobsid = NiObsID(args.obsid_path,args.fyaml_param)
	niobsid.run(args.flag_recreate)



