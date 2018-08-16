#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 August 16'
__version__ = '1.01'
__keyword__ = 'niqlcheck_v%s' % __version__

import os 
import sys 
import glob 
import yaml 
import argparse
import subprocess
import pandas as pd 
import astropy.io.fits as pyfits 

NUM_OF_FPM_ALL = 52 

def bin_spec(src_pha,bgd_pha,
	min_significance=5,max_bins=120,
	emin=0.2,emax=15.0):
	basename = '%s_grp' % os.path.splitext(os.path.basename(src_pha))[0]
	grp_qdp = '%s.qdp' % basename
	grp_pha = '%s.pha' % basename
	cmd  = 'xspec <<EOF\n'
	cmd += 'data 1 %s\n' % src_pha
	cmd += 'back 1 %s\n' % bgd_pha
	cmd += 'setplot energy\n'
	cmd += 'ignore 1:**-%.2f %.2f-**\n' % (emin,emax)
	cmd += 'setplot rebin %d %d\n' % (min_significance,max_bins)
	cmd += 'setplot channel\n'
	cmd += 'cpd /xw\n'
	cmd += 'iplot data\n'
	cmd += 'we %s\n' % basename
	cmd += 'exit\n'
	cmd += 'exit\n'
	print(cmd);os.system(cmd)

	cmd  = 'fgrppha.py %s %s' % (src_pha,grp_qdp)
	print(cmd);os.system(cmd)

	cmd  = 'mv %s %s;' % (grp_pha,os.path.dirname(src_pha))
	cmd += 'rm -f %s.{qdp,pco};' % basename
	print(cmd);os.system(cmd)

	grp_pha = '%s/%s' % (os.path.dirname(src_pha),grp_pha)
	return grp_pha

class NicerQLData():
	def __init__(self,qldata_path,fname_yaml):	
		self.qldata_path = qldata_path
		print("qldata_path: %s" % self.qldata_path)
		
		self.fname_yaml = fname_yaml
		if not os.path.exists(self.fname_yaml):
			print("error: input yaml file does not exists, %s." % self.fname_yaml)			
			exit()
		print("fname_yaml: %s" % self.fname_yaml)
		self.param = yaml.load(open(self.fname_yaml))

	def check_files(self):
		for file in glob.glob('%s/*' % self.qldata_path):
			if os.path.splitext(file)[-1] == '.mkf':
				self.mkf = file
				print("mkf: %s" % self.mkf)
			if os.path.splitext(file)[-1] == '.orb':
				self.orb = file
				print("orb: %s" % self.orb)				
			if os.path.splitext(file)[-1] == '.evt':
				if 'mpu7cl' in file:
					self.clevt = file
					print("clevt: %s" % self.clevt)

					hdu = pyfits.open(self.clevt)
					self.qlobsid = hdu[0].header['OBS_ID']
					self.target  = hdu[0].header['OBJECT']
					self.dateobs = hdu[0].header['DATE-OBS']
					self.exposure = hdu[1].header['EXPOSURE']
					print("qlobsid: %s" % self.qlobsid)

	def make_directory(self):
		self.outdir = '%s/%s' % (self.param['outdir_root_name'],self.qlobsid)
		cmd = 'rm -rf %s; mkdir -p %s' % (self.outdir,self.outdir)
		print(cmd);os.system(cmd)

	def exclude_noisy_fpm(self):
		num_of_fpm = NUM_OF_FPM_ALL - len(self.param['exclude_detid_list'])
		self.cl2evt = '%s/%s_nfpm%d.evt' % (self.outdir,os.path.splitext(os.path.basename(self.clevt))[0],num_of_fpm)
		cmd  = 'fselect %s %s ' % (self.clevt,self.cl2evt)
		cmd += '"'
		for detid in self.param['exclude_detid_list']:
			cmd += '(DET_ID != %s)' % (detid)
			if self.param['exclude_detid_list'].index(detid) < len(self.param['exclude_detid_list'])-1:
				cmd += ' && '
		cmd += '"'		
		print(cmd);os.system(cmd)

	def extract_average_spectrum(self):
		self.cl2pha = self.cl2evt.replace('.evt','.pha')
		cmd  = 'fxselect_extract_spec.py '
		cmd += '-i %s ' % self.cl2evt
		cmd += '-o %s ' % self.cl2pha
		cmd += '-r %s ' % self.param['rmffile']
		cmd += '-a %s ' % self.param['arffile']
		print(cmd);os.system(cmd)

	def fit_average_spectrum(self):
		self.bgd_pha = '%s/niqlcheck/%s' % (os.getenv('NICER_SOFT_PATH'),self.param['bgdfile'])
		self.src_grp_pha = bin_spec(self.cl2pha,self.bgd_pha,
			min_significance=self.param['binspec_min_significance'],
			max_bins=self.param['binspec_max_bins'],
			emin=self.param['binspec_emin'],
			emax=self.param['binspec_emax'])

		fit_model_xcm = '%s/niqlcheck/%s' % (os.getenv('NICER_SOFT_PATH'),self.param['fitmodel'])
		self.fitxcm = self.cl2evt.replace('.evt','_fit.xcm')
		self.fitps = self.cl2evt.replace('.evt','_fit.ps')		
		self.title = '%s %s %s (%.1f s)' % (self.target,self.qlobsid,self.dateobs, self.exposure)
		cmd  = 'xspec <<EOF\n'
		cmd += 'data 1 %s\n' % self.src_grp_pha
		cmd += 'back 1 %s\n' % self.bgd_pha
		cmd += 'setplot energy\n'
		cmd += 'ignore 1:**-%.2f %.2f-**\n' % (self.param['fitspec_emin'],self.param['fitspec_emax'])
		cmd += '@%s\n' % fit_model_xcm
		cmd += 'renorm\n'
		cmd += 'query yes\n'
		cmd += 'fit\n'
		cmd += 'save all %s\n' % self.fitxcm
		cmd += 'cpd /xw\n'
		cmd += 'iplot ld del\n'
		cmd += 'lwid 5\n'
		cmd += 'lwid 5 on 1..100\n'
		cmd += 'csize 1.2\n'
		cmd += 'lab pos y 2.3\n'
		cmd += 'la t %s\n' % self.title
		cmd += 'hard %s/cps\n' % self.fitps 
		cmd += 'exit\n'
		cmd += 'exit\n'
		print(cmd);os.system(cmd)		

		self.fitpdf = os.path.basename(self.fitps).replace('.ps','.pdf')
		cmd  = 'ps2pdf %s;' % self.fitps 
		cmd += 'mv %s %s;' % (self.fitpdf,self.outdir)
		cmd += 'rm -f %s' % self.fitps
		print(cmd);os.system(cmd)

	def select_energy_band(self):
		self.enesel_event_list = []
		for energy_band in self.param['energy_range_list']:
			pi_min = int(float(energy_band[0]) * KEV_TO_PI)
			pi_max = int(float(energy_band[1]) * KEV_TO_PI)
			print(pi_min,pi_max)
			add_str_ext = '_%sto%skeV.evt' % (str(energy_band[0]).replace('.','p'),str(energy_band[1]).replace('.','p'))
			outevt = self.cl2evt.replace('.evt',add_str_ext)
			self.enesel_event_list.append(outevt)
			cmd  = 'fselect %s %s ' % (self.cl2evt,outevt)
			cmd += '"( PI >= %d ) && ( PI < %d )"' % (pi_min,pi_max)
			print(cmd);os.system(cmd)

	def plot_light_curve(self,select_band=[0,1]):
		if len(select_band) != 2:
			print("error: invalid select_eband.")
			exit()
		add_str_ext = '_eband%d_%d.flc' % (select_band[0],select_band[1])
		self.lc_hardness = self.cl2evt.replace('.evt',add_str_ext)
		self.lc_hardness_ps = self.lc_hardness.replace('.flc','_flc.ps')
		cmd  = 'lcurve 2 %s %s ' % (self.enesel_event_list[0],self.enesel_event_list[1])
		cmd += 'window="-" dtnb=%d nbint=%d ' % (self.param['lc_dtnb'],self.param['lc_nbint'])
		cmd += 'outfile=%s ' % self.lc_hardness
		cmd += 'plot=yes plotdev=/xw plotdnum=3 <<EOF\n'
		cmd += 'lab G2 cps (%.1f-%.1f keV)\n' % (self.param['energy_range_list'][select_band[0]][0],self.param['energy_range_list'][select_band[0]][1])
		cmd += 'lab G3 cps (%.1f-%.1f keV)\n' % (self.param['energy_range_list'][select_band[1]][0],self.param['energy_range_list'][select_band[1]][1])		
		cmd += 'lab G4 Hardness \n'
		cmd += 'col 2 on 2\n'
		cmd += 'col 3 on 3\n'		
		cmd += 'col 4 on 4\n'		
		cmd += 'lwid 5\n'		
		cmd += 'hard %s/cps\n' % self.lc_hardness_ps
		cmd += 'exit\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)

		cmd = 'ps2pdf %s' % self.lc_hardness_ps
		print(cmd);os.system(cmd)
		self.lc_hardness_pdf = os.path.basename(self.lc_hardness_ps).replace('.ps','.pdf')
		cmd = 'mv %s %s' % (self.lc_hardness_pdf,self.outdir)
		print(cmd);os.system(cmd)

	def correct_gti(self):
		hdu = pyfits.open(self.cl2evt)
		original_gti_data = hdu['GTI'].data
		num_of_original_gti = len(original_gti_data)
		self.corrected_gti_list = []
		for i in range(num_of_original_gti):
			tstart_i = float(original_gti_data[i][0])
			tstop_i  = float(original_gti_data[i][1])
			if i == 0:
				self.corrected_gti_list.append([tstart_i,tstop_i])
				continue 
			if tstart_i < original_gti_data[i-1][1] + self.param['gti_corr_threshold']:
				self.corrected_gti_list[-1][1] = tstop_i
			else:
				self.corrected_gti_list.append([tstart_i,tstop_i])

		for gti in self.corrected_gti_list:
			print(gti)

	def devide_to_gti(self):
		self.csv_data_list = []
		for gtinum in range(len(self.corrected_gti_list)):		
			tstart = self.corrected_gti_list[gtinum][0]
			tstop = self.corrected_gti_list[gtinum][1]

			subdir = '%s/gti/gti%d' % (self.outdir,gtinum)
			cmd = 'rm -rf %s; mkdir -p %s' % (subdir,subdir)
			print(cmd);os.system(cmd)

			basename = os.path.splitext(os.path.basename(self.cl2evt))[0]
			outevt = '%s/%s_gti%d.evt' % (subdir,basename,gtinum)
			cmd  = 'fxselect_filter_time.py -i %s -o %s ' % (self.cl2evt,outevt)
			cmd += '-d %.7f -u %.7f' % (tstart,tstop)
			print(cmd);os.system(cmd)		

			outpha = outevt.replace('.evt','.pha')
			cmd  = 'fxselect_extract_spec.py '
			cmd += '-i %s ' % outevt
			cmd += '-o %s ' % outpha
			cmd += '-r %s ' % self.param['rmffile']
			cmd += '-a %s ' % self.param['arffile']
			print(cmd);os.system(cmd)

			bgdpha = ""
			data_id = "gti%d" % gtinum
			self.csv_data_list.append([data_id, outpha,bgdpha,self.param['rmffile'],self.param['arffile']])

		self.fcsvfile = '%s/%s_gti.csv' % (self.outdir,self.qlobsid)
		df = pd.DataFrame(self.csv_data_list,columns=['data_id','src_pha','bgd_pha','rmffile','arffile'])
		df.to_csv(self.fcsvfile)
		
	def run(self):
		self.check_files()
		self.make_directory()
		self.exclude_noisy_fpm()
		self.extract_average_spectrum()
		self.fit_average_spectrum()
		#self.correct_gti()
		#self.devide_to_gti()
		#self.select_energy_band()
		#self.plot_light_curve()

if __name__=="__main__":

	sys.stdout.write('... run a single ObsID QL check (%s) ...\n' % __keyword__)

	cmd = 'which fhelp'
	resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if len(resp.stdout.readlines()) == 0:
		sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
		exit()

	parser = argparse.ArgumentParser(
		prog='niqlcheck.py',
		usage='niqlcheck.py obsid_path [--fyaml_param]',
		description='A script for NICER QL data check. [output=spec]',
		epilog='',
		add_help=True,
		)

	NICER_DEFAULT_QLPARAM_FILE = '%s/niqlcheck/niqlcheck_setup_default.yaml' % os.getenv('NICER_SOFT_PATH')

	parser.add_argument(
		'obsid_path',metavar='obsid_path',type=str,        
		help='Input ObsID directory path. The directory should be a single NICER observation directory, which in turn contains xti/{events_uf,events_cl,hk,auxil} subdirectories.')
	parser.add_argument(
		'--fyaml_param',metavar='fyaml_param',type=str,     
		default=NICER_DEFAULT_QLPARAM_FILE,   
		help='yaml file for input parameters.')
	args = parser.parse_args()	
	
	niqldata = NicerQLData(args.obsid_path,args.fyaml_param)
	niqldata.run()
