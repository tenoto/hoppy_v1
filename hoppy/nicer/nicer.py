# -*- coding: utf-8 -*-

import os 
import sys
import glob
import yaml 
import pandas as pd
import inspect
import numpy as np
from time import sleep

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.time import Time
from astropy.stats import bayesian_blocks

def nibackgen3C50(outdir,basename,totspec_name,bkgspec_name,
		bkgidxdir,bkglibdir,gainepoch,
		mode='obsid',obsid=None,clfile=None,ufafile=None,
		dtmin=20.0,dtmax=60.6,hbgcut=0.5,rmf=None,arf=None):
	print("-----[func:%s]" % sys._getframe().f_code.co_name)

	suboutdir = '%s/spec' % outdir
	cmd  = 'rm -rf %s;\n' % suboutdir
	cmd += 'mkdir -p %s\n' % suboutdir
	print(cmd);os.system(cmd)

	flog = '%s/nibackgen3C50_%s.log' % (suboutdir,basename)

	if mode=='obsid':
		clfile = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (outdir,obsid)
	hdu = fits.open(clfile)
	if len(hdu['EVENTS'].data)==0:
		f = open(flog,'w')
		message = 'no cleaned event in the original file %s' % self.clevt
		print(message)
		f.write(message)
		f.close()
		return -1

	cmd = 'rm -f %s.pi %s.pi;' % (totspec_name,bkgspec_name)
	print(cmd);os.system(cmd)

	# prepare a script for nicerl2 
	fcmd = '%s/nibackgen3C50_%s.sh' % (suboutdir,basename)
	f = open(fcmd,'w')
	dump  = '#!/bin/sh -f\n'
	dump += 'nibackgen3C50 '
	dump += 'totspec=\'%s\' ' % totspec_name
	dump += 'bkgspec=\'%s\' ' % bkgspec_name	
	dump += 'bkgidxdir=\'%s\' ' % bkgidxdir
	dump += 'bkglibdir=\'%s\' ' % bkgidxdir	
	dump += 'gainepoch=\'%s\' ' % gainepoch
	dump += 'dtmin=%.1f ' % dtmin
	dump += 'dtmax=%.1f ' % dtmax
	dump += 'hbgcut=%.1f ' % hbgcut
	if mode == 'obsid':
		dump += 'rootdir=\'%s\' ' % outdir.replace(obsid,'')
		dump += 'obsid=\'%s\' ' % obsid
	elif mode == 'segment':
		dump += 'rootdir=\'NONE\' '
		dump += 'obsid=\'NONE\' ' 
		dump += 'calevtdir=\'NONE\' ' 		
		dump += 'clfile=\'%s\' ' % clfile
		dump += 'ufafile=\'%s\' ' % ufafile 
	dump += 'clobber=yes '
	dump += 'chatter=5 '
	dump += '> %s 2>&1 ' % flog
	dump += '\n'
	f.write(dump)
	f.close()	

	# run the script.
	cmd  = 'chmod +x %s\n' % fcmd
	cmd += './%s' % fcmd
	print(cmd);os.system(cmd)

	#f = open(flog,'a')
	#f.write(fversion+'\n')
	#f.close()		


	if os.path.exists(totspec_name+'.pi'):
		flag_exists_totspec = True
	else:
		flag_exists_totspec = False
		print("spectral file %s was not generated." % totspec_name)
		return -1

	if os.path.exists(bkgspec_name+'.pi'):
		flag_exists_bkgspec = True
	else:
		flag_exists_bkgspec = False
		print("spectral file %s was not generated." % bkgspec_name)
		return -1			

	cmd = 'mv %s*.pi %s;' % (totspec_name,suboutdir)
	cmd += 'mv %s*.pi %s;' % (bkgspec_name,suboutdir)	
	print(cmd);os.system(cmd)		

	totspec = '%s/%s.pi' % (suboutdir,totspec_name)
	bkgspec = '%s/%s.pi' % (suboutdir,bkgspec_name)

	if rmf != None:
		cmd  = 'fparkey %s %s RESPFILE\n' % (rmf,totspec)
		cmd += 'fparkey %s %s RESPFILE\n' % (rmf,bkgspec)
		print(cmd);os.system(cmd)
	if arf != None:
		cmd  = 'fparkey %s %s ANCRFILE\n' % (arf,totspec)
		cmd += 'fparkey %s %s ANCRFILE\n' % (arf,bkgspec)
		print(cmd);os.system(cmd)	

	print("-----[func:%s] DONE" % sys._getframe().f_code.co_name)
	return 0

def plot_xspec(totspec,bkgspec,
	xspec_rebin_sigma=5,xspec_rebin_maxnum=100,
	xspec_emin=0.2,xspec_emax=10.0,xspec_plot_type='ld',
	xspec_ymin=1e-4,xspec_ymax=100):
	print("-----[func:%s]" % sys._getframe().f_code.co_name)

	print(totspec)
	print(bkgspec)	
	hdu = fits.open(totspec)
	keyword_dateobs = hdu[0].header['DATE-OBS']
	keyword_object = hdu[0].header['OBJECT']		
	keyword_obsid = hdu[0].header['OBS_ID']
	keyword_exposure = hdu[0].header['EXPOSURE']			
	title = '%s ObsID:%s (%s, %.1f sec)' % (keyword_object, keyword_obsid, keyword_dateobs,keyword_exposure)
	fspec = os.path.splitext(os.path.basename(totspec))[0]
	cmd  = "xspec <<EOF\n"
	cmd += "data 1 %s\n" % totspec
	cmd += "data 2 %s\n" % totspec
	cmd += "back 2 %s\n" % bkgspec
	cmd += "data 3 %s\n" % bkgspec	
	cmd += "back 3 none\n"
	cmd += "setplot energy\n"
	cmd += "setplot rebin %d %d\n" % (xspec_rebin_sigma,xspec_rebin_maxnum)
	cmd += "ignore 1-3:**-%.1f,%.1f-**\n" % (xspec_emin,xspec_emax)
	cmd += "iplot %s\n" % xspec_plot_type
	cmd += "time off\n"
	cmd += "lwid 5\n"
	cmd += "lwid 5 on 1..3\n"
	cmd += "la t %s\n" % title
	cmd += "r x %.1f %.1f\n" % (xspec_emin,xspec_emax)
	#cmd += "r y 1e-4 100.0\n"
	cmd += "r y %s %s\n" % (xspec_ymin,xspec_ymax)
	cmd += "hard %s.ps/cps\n" % fspec
	cmd += "exit\n"
	cmd += "exit\n"
	cmd += "EOF\n"
	print(cmd)
	os.system(cmd)

	outdir = os.path.dirname(totspec)
	cmd  = "ps2pdf %s.ps\n" % fspec
	cmd += "rm -f %s.ps\n" % fspec
	cmd += "mv %s.pdf %s/\n" % (fspec,outdir)
	print(cmd)
	os.system(cmd)	

	print("-----[func:%s] Done" % sys._getframe().f_code.co_name)

def plot_lcurve(eventfile,outdir,lc_energy_bands,lc_time_bin_sec):
	print("-----[func:%s]" % sys._getframe().f_code.co_name)

	suboutdir = '%s/lc' % outdir
	cmd  = 'rm -rf %s;\n' % suboutdir
	cmd += 'mkdir -p %s\n' % suboutdir
	print(cmd);os.system(cmd)

	basename = os.path.splitext(os.path.basename(eventfile))[0]
	flog = '%s/lcurve_%s.log' % (suboutdir,basename)
	hdu = fits.open(eventfile)
	if len(hdu['EVENTS'].data)==0:
		f = open(flog,'w')
		message = 'no cleaned event in the original file %s' % self.clevt
		print(message)
		f.write(message)
		f.close()
		return -1

	for eband in lc_energy_bands:
		emin,emax = lc_energy_bands[eband]
		fenesel = '%s_%s.evt' % (basename,eband)
		cmd += "fselect_energy.py %s %s %.1f %.1f;\n" % (eventfile,fenesel,emin,emax)
		cmd += "mv %s %s;\n" % (fenesel,suboutdir)
	print(cmd); os.system(cmd)


	tstart_met = hdu['EVENTS'].data['TIME'][0]
	tstop_met = hdu['EVENTS'].data['TIME'][-1]	
	onsource_sec = tstop_met - tstart_met	
	nbint = round(float(onsource_sec) / float(lc_time_bin_sec))

	keyword_dateobs = hdu[0].header['DATE-OBS']
	keyword_object = hdu[0].header['OBJECT']		
	keyword_obsid = hdu[0].header['OBS_ID']
	title = '%s ObsID:%s (%s, %.1f sec)' % (keyword_object, keyword_obsid, keyword_dateobs, onsource_sec)

	flcfile = '%s/%s_ene.flc' % (suboutdir,basename)
	fcmd = '%s/%s_ene.sh' % (suboutdir,basename)
	#flog = '%s/%s_ene.log' % (suboutdir,basename)
	print(fcmd)
	f = open(fcmd,'w')
	dump  = '#!/bin/sh -f\n'
	dump += 'lcurve nser=%d ' % len(lc_energy_bands)
	i = 1
	for eband in lc_energy_bands:
		fenesel = '%s_%s.evt' % (basename,eband)
		dump += 'cfile%d="%s/%s" ' % (i,suboutdir,fenesel)
		i += 1 
	dump += 'window="-" ' 
	dump += 'dtnb=%d ' % lc_time_bin_sec
	dump += 'nbint=%d ' % nbint
	dump += 'outfile="%s" ' % flcfile 
	dump += 'plotdnum=%s ' % len(lc_energy_bands)
	dump += 'plot=yes plotdev="/null" <<EOF\n'
	dump += 'lwid 5\n'
	#dump += 'la ot %s\n' % title 
	dump += 'lab rotate\n'
	dump += 'lab pos y 2.8\n'
	i = 2
	for eband in lc_energy_bands:
		emin,emax = lc_energy_bands[eband]		
		dump += 'lab y%d %.1f-%.1f keV\n' % (i,emin,emax)
		dump += 'col %d on %d\n' % (i,i)
		i += 1 
	dump += 'hard %s_ene.ps/cps\n' % basename
	dump += 'quit\n'
	dump += 'EOF'
	print(dump)
	f.write(dump)
	f.close()

	# run the script.
	cmd  = 'chmod +x %s\n' % fcmd
	cmd += './%s' % fcmd
	print(cmd);	os.system(cmd)

	cmd  = "ps2pdf %s_ene.ps\n" % basename
	cmd += "rm -f %s_ene.ps\n" % basename
	cmd += "mv %s_ene.pdf %s\n" % (basename,suboutdir)
	print(cmd);	os.system(cmd)		

	print("-----[func:%s] Done" % sys._getframe().f_code.co_name)

def barycen(evtfile,outfile,outdir,orbfile,ra,dec,overwrite=False):
	print("-----[func:%s]" % sys._getframe().f_code.co_name)

	if overwrite:
		cmd  = 'rm -rf %s;\n' % outdir
		cmd += 'mkdir -p %s\n' % outdir
	else:
		cmd = 'mkdir -p %s\n' % outdir
	print(cmd);os.system(cmd)

	baryevt = '%s/%s' % (outdir,outfile)
	basename = os.path.splitext(os.path.basename(outfile))[0]
	fcmd = '%s/barycorr_%s.sh' % (outdir,basename)
	flog = '%s/barycorr_%s.log' % (outdir,basename)
	f = open(fcmd,'w')
	dump  = '#!/bin/sh -f\n'
	dump += 'barycorr infile=%s ' % evtfile
	dump += 'outfile=%s ' % baryevt
	dump += 'orbitfiles=%s ' % orbfile 
	dump += 'ra=%.6f dec=%.6f ' % (ra,dec)
	dump += 'refframe=ICRS ephem=JPLEPH.430 '
	dump += '> %s 2>&1 ' % flog
	dump += '\n'
	f.write(dump)
	f.close()	

	# run the script.
	cmd  = 'chmod +x %s\n' % fcmd
	cmd += './%s' % fcmd
	print(cmd);os.system(cmd)		
	print("-----[func:%s] Done" % sys._getframe().f_code.co_name)
	return 0 

def fit_xspec(phafile,indir,backgrnd,rmffile,arffile,modelxcm,
	binminsig,binmaxbin,fitemin,fitemax,ratebands,fluxbands,parerrnum):
	print("-----[func:%s]" % sys._getframe().f_code.co_name)

	outdir = '%s/fit' % indir
	cmd = 'rm -rf %s' % outdir
	print(cmd);os.system(cmd)

	cmd  = 'xspec_fit.py '
	cmd += '%s ' % phafile
	cmd += '-o %s ' % outdir
	cmd += '-b %s ' % backgrnd
	cmd += '-r %s ' % rmffile
	cmd += '-a %s ' % arffile
	cmd += '-m %s ' % modelxcm
	cmd += '-s %d -n %d ' % (binminsig,binmaxbin)
	cmd += '--fitemin %.1f --fitemax %.1f ' % (fitemin,fitemax)
	cmd += '--rateband %s ' % ratebands
	cmd += '--fluxband %s ' % fluxbands
	cmd += '--parerrnum %s \n' % parerrnum 
	print(cmd);os.system(cmd)

	basename = os.path.splitext(phafile)[0]
	fit_pdf = '%s_s%dn%d_fit.pdf' % (basename,binminsig,binmaxbin)
	fit_pdf = fit_pdf.replace('/spec/','/fit/')

	print("-----[func:%s] Done" % sys._getframe().f_code.co_name)
	return fit_pdf	

class NicerObsID():
	def __init__(self,obsid,indir,outdir,param):
		self.param = param; 
		self.obsid = obsid; 
		self.indir = indir; 
		self.outdir = outdir; 
		print("**[NicerObsID %s] A new NicerObs is generated" % self.obsid)

		self.clevt = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (self.outdir,self.obsid)		
		self.ufaevt = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt' % (self.outdir,self.obsid)		
		self.orbitfile = "%s/auxil/ni%s.orb.gz" % (self.outdir,self.obsid)
		self.totspec = '%s/spec/ni%s_3c50_ave_tot.pi' % (self.outdir,self.obsid)
		self.bkgspec = '%s/spec/ni%s_3c50_ave_bkg.pi' % (self.outdir,self.obsid)		

		self.total_number_of_block = 0

		self.rate = None
		self.rate_error = None
		self.band_low = None
		self.band_high = None

		self.setup_yamlfile = '%s/ni%s_result.yaml' % (self.outdir,self.obsid)
		self.have_event_files()
		self.set_basic_information()
		sleep(1) # this is needed, otherwise sometimes flag_clevt_has_events can not be ste.
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))

	def set_basic_information(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))	
		if os.path.exists(self.clevt):
			hdu = fits.open(self.clevt)
			self.date_obs = hdu['EVENTS'].header['DATE-OBS']
			self.object = hdu['EVENTS'].header['OBJECT']
			self.exposure = float(hdu['EVENTS'].header['EXPOSURE'])
		else:
			self.date_obs = None
			self.object = None
			self.exposure = None
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0 

	def set_nisegment_list(self):	
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))	
		print(hasattr(self, 'nigsegment_list'))
		if not hasattr(self, 'nigsegment_list') and os.path.exists('%s/segment' % self.outdir):
			self.dir_segment_main = '%s/segment' % self.outdir
			print(self.dir_segment_main)
			self.number_of_segment = len(glob.glob('%s/seg*' % self.dir_segment_main))
			print(self.number_of_segment)
			self.nigsegment_list = []
			for i in range(self.number_of_segment):
				segment_num = i + 1 
				segment_dir = '%s/seg%003d' % (self.dir_segment_main,segment_num)			
				print("----- segment_dir = %s" % segment_dir)			
				segment_basename = 'ni%s_0mpu7_seg%03d' % (self.obsid,segment_num)
				print(segment_num)			
				print(segment_dir)
				print(segment_basename)

				gtifile = '%s/ni%s_0mpu7_seg%03d.gti' % (segment_dir,self.obsid,segment_num)
				hdu_gti = fits.open(gtifile)
				tstart  = hdu_gti['GTI'].data[0]['START']
				tstop  = hdu_gti['GTI'].data[0]['STOP']
				print(tstart,tstop)
				self.nigsegment_list.append(NicerInterval(
					dtype="segment",
					parent_nicerobsid=self,
					parent_segment=None,
					num=segment_num,
					outdir=segment_dir,
					basename=segment_basename,
					tstart=tstart,
					tstop=tstop))
		print(self.nigsegment_list)
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0 

	def have_event_files(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		self.flag_clevt_has_events = False
		if os.path.exists(self.clevt):
			self.flag_exists_clevt = True
			hdu = fits.open(self.clevt)
			self.num_of_clevt = len(hdu['EVENTS'].data)			
			if self.num_of_clevt > 0:
				self.flag_clevt_has_events = True
			else:
				self.flag_clevt_has_events = False
		else:
			self.flag_exists_clevt = False
			self.num_of_clevt = None

		if os.path.exists(self.ufaevt):
			self.flag_exists_ufaevt = True
		else:
			self.flag_exists_ufaevt = False
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0 

	def make_directory(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		cmd  = 'rm -rf %s;\n' % (self.outdir)
		cmd += 'mkdir -p %s/xti;\n' % (self.outdir)
		cmd += 'ln -s %s/{auxil,log} %s/;\n' % (self.indir,self.outdir)
		cmd += 'ln -s %s/xti/{event_uf,hk} %s/xti/;\n' % (self.indir,self.outdir)	
		print(cmd);os.system(cmd)
		#self.dump_setup_to_yamlfile()
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		#return 0 

	def dump_setup_to_yamlfile(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		print("-----------")
		output_dict = {}
		for attr in dir(self):
			if not attr.startswith('__') and not inspect.ismethod(getattr(self,attr)):
				print("self.%s = %r" % (attr, getattr(self, attr))) 
				output_dict[attr] = getattr(self, attr)
		print("-----------")			
		cmd = 'rm -f %s;\n' % self.setup_yamlfile
		print(cmd);os.system(cmd)
		with open(self.setup_yamlfile, "a") as wf:
			yaml.dump(output_dict, wf, default_flow_style=False)		
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))

	def run_nicerl2(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		fversion = os.popen('fversion').read().rstrip()		
	
		# prepare a script for nicerl2 for each ObsID
		self.fcmd_nicerl2 = '%s/nicerl2_%s.sh' % (self.outdir,self.obsid)
		self.flog_nicerl2 = '%s/nicerl2_%s.log' % (self.outdir,self.obsid)
		f = open(self.fcmd_nicerl2,'w')
		dump  = '#!/bin/sh -f\n'
		dump += 'nicerl2 indir=%s ' % (self.outdir)
		dump += 'picalfile=%s ' % self.param['nicerl2_gcalfile']
		dump += 'clobber=yes '
		dump += '> %s 2>&1 ' % self.flog_nicerl2
		dump += '\n'
		f.write(dump)
		f.close()

		# run the script.
		cmd  = 'chmod +x %s\n' % self.fcmd_nicerl2
		cmd += './%s' % self.fcmd_nicerl2
		print(cmd)
		os.system(cmd)

		f = open(self.flog_nicerl2,'a')
		f.write(fversion+'\n')
		f.close()

		self.have_event_files()
		self.set_basic_information()		
		self.dump_setup_to_yamlfile()	
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def run_nibackgen3C50(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		fversion = os.popen('fversion').read().rstrip()		

		#if not self.flag_clevt_has_events:
		#	print("No cleaned events.")
		#	return 0

		nibackgen3C50(
			outdir=self.outdir,
			basename='%s' % self.obsid,
			totspec_name='ni%s_3c50_ave_tot' % self.obsid,
			bkgspec_name='ni%s_3c50_ave_bkg' % self.obsid,
			bkgidxdir=self.param['nibackgen3c50_bkgidxdir'],
			bkglibdir=self.param['nibackgen3c50_bkgidxdir'],
			gainepoch=self.param['nibackgen3c50_gainepoch'],
			mode='obsid',
			obsid=self.obsid,
			dtmin=self.param['nibackgen3c50_dtmin'],
			dtmax=self.param['nibackgen3c50_dtmax'],
			hbgcut=self.param['nibackgen3c50_hbgcut'],
			rmf=self.param['xspec_rmf'],
			arf=self.param['xspec_arf'])

		plot_xspec(
			totspec=self.totspec,
			bkgspec=self.bkgspec,
			xspec_rebin_sigma=self.param['xspec_rebin_sigma'],
			xspec_rebin_maxnum=self.param['xspec_rebin_maxnum'],
			xspec_emin=self.param['xspec_emin'],
			xspec_emax=self.param['xspec_emax'],
			xspec_plot_type=self.param['xspec_plot_type'],
			xspec_ymin=self.param['xspec_ymin'],
			xspec_ymax=self.param['xspec_ymax'])
		self.spec_pdf = '%s/spec/ni%s_3c50_ave_tot.pdf' % (self.outdir,self.obsid)

		self.dump_setup_to_yamlfile()
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def fit_avespec(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		fversion = os.popen('fversion').read().rstrip()		

		#if not self.flag_clevt_has_events:
		#	print("No cleaned events.")
		#	return 0

		self.fit_avespec_pdf = fit_xspec(
			phafile=self.totspec,
			indir=self.outdir,
			backgrnd=self.bkgspec,
			rmffile=self.param['xspec_rmf'],
			arffile=self.param['xspec_arf'],
			modelxcm=self.param['xspec_model'],
			binminsig=self.param['xspec_rebin_sigma'],binmaxbin=self.param['xspec_rebin_maxnum'],
			fitemin=self.param['xspec_emin'],fitemax=self.param['xspec_emax'],
			ratebands=self.param['xspec_rateband'],fluxbands=self.param['xspec_fluxband'],
			parerrnum=str(self.param['xspec_parerrnum']).replace(' ',''))
		print(self.fit_avespec_pdf)

		fit_yaml = self.fit_avespec_pdf.replace(".pdf",".yaml")
		if not os.path.exists(fit_yaml):
			print("warning: %s does not exist." % fit_yaml)
			return 0
		fit_param = yaml.load(open(fit_yaml),Loader=yaml.FullLoader)
		self.rate = fit_param['ratelist'][0][0]
		self.rate_error = fit_param['ratelist'][0][1]
		self.band_low = fit_param['ratebands'][0][0] 
		self.band_high = fit_param['ratebands'][0][1] 		

		self.dump_setup_to_yamlfile()
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def plot_lightcurve(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))


		#if not self.flag_clevt_has_events:
		#	print("No cleaned events.")
		#	return 0

		plot_lcurve(
			eventfile=self.clevt,
			outdir=self.outdir,
			lc_energy_bands=self.param['lc_energy_bands'],
			lc_time_bin_sec=self.param['lc_time_bin_sec']
			)

		self.lc_pdf = '%s/lc/ni%s_0mpu7_cl_ene.pdf' % (self.outdir,self.obsid)

		self.dump_setup_to_yamlfile()
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0		

	def devide_to_segment(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))


		if not self.flag_clevt_has_events:
			print("No cleaned events.")
			return 0

		self.dir_segment_main = '%s/segment' % self.outdir
		cmd  = 'rm -rf %s;\n' % self.dir_segment_main
		cmd += 'mkdir -p %s\n' % self.dir_segment_main
		print(cmd);os.system(cmd)

		# re-arrange GTIs
		cmd = 'rm -f tmp_gti.txt; group_gti.py %s tmp_gti.txt --gti_separation_threshold %.2f' % (self.clevt,self.param['gti_separation_threshold'])
		print(cmd);os.system(cmd)

		# convert text file GTI format to fits file.
		segment_gti = '%s/ni%s_segment.gti' % (self.dir_segment_main,self.obsid)
		cmd = 'fconv_txt2gti.py -i tmp_gti.txt -o %s; rm -f tmp_gti.txt' % segment_gti 
		print(cmd);os.system(cmd)

		hdu = fits.open(segment_gti)
		self.number_of_segment = len(hdu['GTI'].data)
		print("Number of segments = %d" % self.number_of_segment)

		self.nigsegment_list = []
		for i in range(self.number_of_segment):
			segment_num = i + 1 
			segment_dir = '%s/seg%03d' % (self.dir_segment_main,segment_num)
			segment_basename = 'ni%s_0mpu7_seg%03d' % (self.obsid,segment_num)
			segment_tstart = hdu['GTI'].data[i][0]
			segment_tstop = hdu['GTI'].data[i][1]
			print(segment_tstart,segment_tstop)
			self.nigsegment_list.append(NicerInterval(
				dtype="segment",
				parent_nicerobsid=self,
				parent_segment=None,
				num=segment_num,
				outdir=segment_dir,
				basename=segment_basename,
				tstart=segment_tstart,
				tstop=segment_tstop))
			self.nigsegment_list[-1].make_directory()
			self.nigsegment_list[-1].filter_event()
			self.nigsegment_list[-1].have_event_files()
			if self.nigsegment_list[-1].flag_exists_clevt and self.nigsegment_list[-1].flag_clevt_has_events:
				self.nigsegment_list[-1].plot_lightcurve()		
				self.nigsegment_list[-1].run_nibackgen3C50()
			else:
				continue

		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def run_barycorr(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))


		#if not self.flag_clevt_has_events:
		#	print("No cleaned events.")
		#	return 0

		basename = os.path.splitext(os.path.basename(self.clevt))[0]
		if self.param['barycorr_RA_deg'] == "None":
			self.param['barycorr_RA_deg'] = float(hdu['EVENTS'].header['RA_OBJ'])
		if self.param['barycorr_DEC_deg'] == "None":
			self.param['barycorr_DEC_deg'] = float(hdu['EVENTS'].header['DEC_OBJ'])		
		barycen(
			evtfile=self.clevt,
			outfile=basename+'_bary.evt',
			outdir=self.outdir+'/bary',
			orbfile=self.orbitfile,
			ra=self.param['barycorr_RA_deg'],dec=self.param['barycorr_DEC_deg'],
			overwrite=True)

		self.dump_setup_to_yamlfile()
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0		

	def fit_of_segment(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		if not self.flag_clevt_has_events:
			print("No cleaned events.")
			return 0

		self.set_nisegment_list()
		print(self.nigsegment_list)
		print("================")
		for nisegment in self.nigsegment_list:
			print(nisegment)

			nisegment.fit_avespec_pdf = fit_xspec(
				phafile=nisegment.totspec_pha,
				indir=nisegment.outdir,
				backgrnd=nisegment.bkgspec_pha,
				rmffile=self.param['xspec_rmf'],
				arffile=self.param['xspec_arf'],
				modelxcm=self.param['xspec_model'],
				binminsig=self.param['xspec_rebin_sigma'],binmaxbin=self.param['xspec_rebin_maxnum'],
				fitemin=self.param['xspec_emin'],fitemax=self.param['xspec_emax'],
				ratebands=self.param['xspec_rateband'],fluxbands=self.param['xspec_fluxband'],
				parerrnum=str(self.param['xspec_parerrnum']).replace(' ',''))
			print(nisegment.fit_avespec_pdf)

		self.dump_setup_to_yamlfile()

		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def devide_to_block(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		#if not self.flag_clevt_has_events:
		#	print("No cleaned events.")
		#	return 0

		self.set_nisegment_list()
		for nisegment in self.nigsegment_list:
			nisegment.devide_to_block()

			self.total_number_of_block += nisegment.number_of_block

		#self.dump_setup_to_yamlfile()
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def fit_of_block(self):
		print("**[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		#if not self.flag_clevt_has_events:
		#	print("No cleaned events.")
		#	return 0

		self.set_nisegment_list()
		for nisegment in self.nigsegment_list:
			for block_dir in glob.glob('%s/block/block*' % nisegment.outdir):
				glob_totspec = glob.glob('%s/spec/*_tot.pi' % (block_dir))
				glob_bkgspec = glob.glob('%s/spec/*_bkg.pi' % (block_dir))
				if not (len(glob_totspec) == 1 and len(glob_bkgspec) == 1):
					print("skip %s" % block_dir)
					continue
				totspec_pha = glob_totspec[0]
				bkgspec_pha = glob_bkgspec[0]				

				fit_avespec_pdf = fit_xspec(
					phafile=totspec_pha,
					indir=block_dir,
					backgrnd=bkgspec_pha,
					rmffile=self.param['xspec_rmf'],
					arffile=self.param['xspec_arf'],
					modelxcm=self.param['xspec_model'],
					binminsig=self.param['xspec_rebin_sigma'],binmaxbin=self.param['xspec_rebin_maxnum'],
					fitemin=self.param['xspec_emin'],fitemax=self.param['xspec_emax'],
					ratebands=self.param['xspec_rateband'],fluxbands=self.param['xspec_fluxband'],
					parerrnum=str(self.param['xspec_parerrnum']).replace(' ',''))
			print(fit_avespec_pdf)

		#self.dump_setup_to_yamlfile()
		print("**[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

class NicerInterval():
	def __init__(self,dtype="segment",parent_nicerobsid=None,parent_segment=None,
		num=None,outdir=None,basename=None,tstart=None,tstop=None):
		self.dtype = dtype
		self.parent_nicerobsid = parent_nicerobsid
		self.parent_segment = parent_segment
		self.num = num
		self.outdir = outdir
		self.basename = basename
		self.tstart = tstart
		self.tstop = tstop

		self.obsid = self.parent_nicerobsid.obsid 
		self.param = self.parent_nicerobsid.param
		if self.parent_segment != None:
			self.segment_num = self.parent_segment.num 
		else:
			self.segment_num = None

		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] A NicerInterval is generated." %(self.dtype,self.obsid,self.segment_num,self.num))

		self.gti_txt = '%s/%s_gti.txt' % (self.outdir,self.basename)
		self.gti = '%s/%s.gti' % (self.outdir,self.basename)
		self.clevt = '%s/%s_cl.evt' % (self.outdir,self.basename)
		self.ufaevt = '%s/%s_ufa.evt' % (self.outdir,self.basename)

		self.totspec_name = '%s_3c50_tot' % (self.basename)
		self.bkgspec_name = '%s_3c50_bkg' % (self.basename)
		self.totspec_pha = '%s/spec/%s.pi' % (self.outdir,self.totspec_name)
		self.bkgspec_pha = '%s/spec/%s.pi' % (self.outdir,self.bkgspec_name)	

		self.spec_pdf = None
		self.lc_pdf = None

		self.duration = self.tstop - self.tstart 
		self.nbint = round(self.duration/float(self.param['lc_time_bin_sec']))
		if dtype == "segment":
			self.title = 'Segment (ObsID=%s Seg=%003d) %.2f-%.2f (%.1f s)' % (self.obsid,self.num,self.tstart,self.tstop,self.duration)
		elif dtype == "block":
			self.title = 'Block (ObsID=%s Seg=%003d Block=%03d) %.2f-%.2f (%.1f s)' % (self.obsid,self.segment_num,self.num,self.tstart,self.tstop,self.duration)

		self.block_dir = None

	def make_directory(self):
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))
		cmd = 'rm -rf %s; mkdir -p %s\n' % (self.outdir,self.outdir)
		print(cmd);os.system(cmd)

		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s DONE" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))
		return 0

	def filter_event(self):
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))

		f = open(self.gti_txt,'w')
		dump = '%10f, %10f' % (self.tstart,self.tstop)
		f.write(dump)
		f.close()
		cmd = 'fconv_txt2gti.py -i %s -o %s' % (self.gti_txt,self.gti)
		print(cmd);os.system(cmd)

		cmd  = 'xselect_gtifilter.py '
		cmd += '-i %s ' % self.parent_nicerobsid.clevt 
		cmd += '-g %s ' % self.gti 
		cmd += '-o %s;\n' % self.clevt
		print(cmd);os.system(cmd)

		cmd  = 'xselect_gtifilter.py '
		cmd += '-i %s ' % self.parent_nicerobsid.ufaevt 
		cmd += '-g %s ' % self.gti 
		cmd += '-o %s;\n' % self.ufaevt
		print(cmd);os.system(cmd)
	
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s DONE" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))
		return 0

	def plot_lightcurve(self):
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))

		plot_lcurve(
			eventfile=self.clevt,
			outdir=self.outdir,
			lc_energy_bands=self.param['lc_energy_bands'],
			lc_time_bin_sec=self.param['lc_time_bin_sec']
			)

		self.lc_pdf = '%s/%s_lc.pdf' % (self.outdir,self.basename)

		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s DONE" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))
		return 0

	def run_nibackgen3C50(self):
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))

		flag_nibackgen3c50 = nibackgen3C50(
			outdir=self.outdir,
			basename=self.basename,
			totspec_name=self.totspec_name,
			bkgspec_name=self.bkgspec_name,
			bkgidxdir=self.param['nibackgen3c50_bkgidxdir'],
			bkglibdir=self.param['nibackgen3c50_bkgidxdir'],
			gainepoch=self.param['nibackgen3c50_gainepoch'],
			mode='segment',
			dtmin=self.param['nibackgen3c50_dtmin'],
			dtmax=self.param['nibackgen3c50_dtmax'],
			hbgcut=self.param['nibackgen3c50_hbgcut'],
			rmf=self.param['xspec_rmf'],
			arf=self.param['xspec_arf'],
			clfile=self.clevt,
			ufafile=self.ufaevt)

		if flag_nibackgen3c50 == 0:
			plot_xspec(
				totspec=self.totspec_pha,
				bkgspec=self.bkgspec_pha,
				xspec_rebin_sigma=self.param['xspec_rebin_sigma'],
				xspec_rebin_maxnum=self.param['xspec_rebin_maxnum'],
				xspec_emin=self.param['xspec_emin'],
				xspec_emax=self.param['xspec_emax'],
				xspec_plot_type=self.param['xspec_plot_type'],
				xspec_ymin=self.param['xspec_ymin'],
				xspec_ymax=self.param['xspec_ymax'])
			self.spec_pdf = '%s/spec/%s.pdf' % (self.outdir,self.totspec_name)

		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s DONE" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))
		return 0

	def have_event_files(self):
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))

		if os.path.exists(self.clevt):
			self.flag_exists_clevt = True
			hdu = fits.open(self.clevt)
			if len(hdu['EVENTS'].data) > 0:
				self.flag_clevt_has_events = True
			else:
				self.flag_clevt_has_events = False
		else:
			self.flag_exists_clevt = False

		if os.path.exists(self.ufaevt):
			self.flag_exists_ufaevt = True
		else:
			self.flag_exists_ufaevt = False
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s DONE" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))
		return 0

	def devide_to_block(self):
		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))

		self.block_dir = '%s/block' % self.outdir
		cmd = 'rm -rf %s; mkdir -p %s;' % (self.block_dir,self.block_dir)
		print(cmd);os.system(cmd)

		self.lcfile = '%s/lc/ni%s_0mpu7_seg%03d_cl_ene.flc' % (self.outdir,self.obsid,self.num)
		if not os.path.exists(self.lcfile):
			print("file does not exist: %s" % self.lcfile)
			return -1

		self.lchdu = fits.open(self.lcfile)
		t = self.lchdu['RATE'].data['TIME']
		x = self.lchdu['RATE'].data['RATE1']
		xe = self.lchdu['RATE'].data['ERROR1']

		filter_x = np.where(~np.isnan(x))
		t = t[filter_x]
		x = x[filter_x]
		xe = xe[filter_x]
		self.edges = bayesian_blocks(t,x,xe,
			fitness=self.parent_nicerobsid.param["bayesian_blocks_fitness"],
			p0=float(self.parent_nicerobsid.param["bayesian_blocks_p0"]))
		print(self.edges)		

		self.timezero_i = float(self.lchdu['RATE'].header['TIMEZERI'])
		self.timezero_f = float(self.lchdu['RATE'].header['TIMEZERF'])
		self.timezero_mjd = self.timezero_i + self.timezero_f + 40000
		print(self.timezero_mjd)

		cmd = 'nitimeconv.py %.10f -f mjd -s utc > tmp.log' % self.timezero_mjd
		print(cmd);os.system(cmd)
		for line in open('tmp.log'):
			cols = line.split()
			if len(cols) == 0:
				continue
			if cols[0] == 'NICER':
				self.met = float(cols[7])
				break
		print(self.met)		
		cmd = 'rm -f tmp.log'
		print(cmd);os.system(cmd)

		self.edges_met = []
		for ed in self.edges:
			self.edges_met.append(ed + self.met)

		self.segment_baysian_blocks_log = '%s/ni%s_seg%03d_edges.log' % (
			self.block_dir,self.obsid,self.num)
		f = open(self.segment_baysian_blocks_log,'w')
		dump = 'NICER ObsID:%s, Segment:%03d\n' % (self.obsid,self.num)
		dump += 'TIMEZERI:%d\n' % self.lchdu['RATE'].header['TIMEZERI']
		dump += 'TIMEZERF:%.10f\n' % self.lchdu['RATE'].header['TIMEZERF']
		dump += 'MJD START:%.10f\n' % self.timezero_mjd
		dump += 'MET START:%.10f\n' % self.met 
		dump += 'Edges: %s\n' % self.edges
		dump += 'Edges (MET): %s\n' % self.edges_met
		f.write(dump)
		f.close()

		self.segment_baysian_blocks_lc_pdf = '%s/ni%s_seg%03d_lc.pdf' % (
			self.block_dir,self.obsid,self.num)
		fig = plt.figure(figsize=(12, 8))
		plt.errorbar(t,x,yerr=xe,marker='o')
		for ed in self.edges:
			plt.axvline(ed)
		plt.savefig(self.segment_baysian_blocks_lc_pdf)

		self.number_of_block = len(self.edges_met) - 1 
		self.niblock_list = []
		for i in range(self.number_of_block):
			block_num = i + 1 

			block_dir = '%s/block%03d' % (self.block_dir,block_num)
			block_basename = 'ni%s_0mpu7_seg%03d_block%03d' % (self.obsid,self.num,block_num)
			block_tstart = self.edges_met[i]
			block_tstop = self.edges_met[i+1]
			print(block_tstart,block_tstop)

			self.niblock_list.append(NicerInterval(
				dtype="block",
				parent_nicerobsid=self.parent_nicerobsid,
				parent_segment=self,
				num=block_num,
				outdir=block_dir,
				basename=block_basename,
				tstart=block_tstart,
				tstop=block_tstop))
			self.niblock_list[-1].make_directory()
			self.niblock_list[-1].filter_event()
			self.niblock_list[-1].have_event_files()
			self.niblock_list[-1].plot_lightcurve()
			self.niblock_list[-1].run_nibackgen3C50()

		print("** [%s ObsID=%s Parent_Seg=%s Num=%s] %s DONE" %(self.dtype,self.obsid,self.segment_num,self.num,sys._getframe().f_code.co_name))
		return 0	

class NicerProcessLog():
	def __init__(self,csvfile):
		print("*[NicerProcessLog] A new NicerProcessLog object is generated.")

		self.csvfile = csvfile
		self.htmlfile = os.path.splitext(self.csvfile)[0] + '.html'

		self.columns = ['object','date_obs','exp(ks)','#evt','dir','spec','rate','fit','lc','seg','#block','mkdir','nicerl2','3c50','xspec','lcurve','bary','div2seg','fit2seg','div2block','fit2block']
		self.dict_null = {}
		for col in self.columns:
			self.dict_null[col] = "Unproc"

		if os.path.exists(self.csvfile):
			print("...reading %s." % self.csvfile )
			self.df = pd.read_csv(csvfile,dtype={'ObsID':str})
			self.df.set_index('ObsID',inplace=True)
		else:
			print("...making a new process log.")
			self.df = pd.DataFrame(
				columns=self.columns+['ObsID'],
				dtype=str
				)
			self.df.set_index('ObsID',inplace=True)
		self.show()

	def show(self):
		print("******* summary ******")
		print(self.df)
		print("**********************")		

	def write_to_csvfile(self):
		print("*[NicerProcessLog] %s" % (sys._getframe().f_code.co_name))
		self.df.to_csv(self.csvfile)		

	def write_to_htmlfile(self):
		print("*[NicerProcessLog] %s" % (sys._getframe().f_code.co_name))
		raw_html = self.df.to_html(open(self.htmlfile,'w'))
		cmd = 'sed -e "s/&lt;/</g" %s > tmp.html' % self.htmlfile
		print(cmd);os.system(cmd)
		cmd = 'sed -e "s/&gt;/>/g" tmp.html > tmp2.html' 
		print(cmd);os.system(cmd)		
		cmd = 'sed -e "s/<td>Done/<td bgcolor="#01FF70">Done/g" tmp2.html > tmp3.html' 
		print(cmd);os.system(cmd)		
		cmd = 'sed -e "s/<td>Error/<td bgcolor="#FF6666">Error/g" tmp3.html > tmp4.html' 
		print(cmd);os.system(cmd)						
		cmd = 'mv tmp4.html %s; rm -f tmp.html tmp?.html' % self.htmlfile
#		cmd = 'mv tmp2.html %s; ' % self.htmlfile
		print(cmd);os.system(cmd)		


class NicerManager():
	def __init__(self,setup_yamlfile,obsid_lstfile):
		print("*[NicerManager] A new NicerManager object is generated.")

		print('...setup_yamlfile: {}'.format(setup_yamlfile))
		print('...obsid_lstfile: {}'.format(obsid_lstfile))
		self.param = yaml.load(open(setup_yamlfile),Loader=yaml.SafeLoader)
		self.obsid_lstfile = obsid_lstfile
		print('...param: %s' % self.param)

		print('...seting nicerobs_lst')
		self.nicerobs_lst = []
		for line in open(self.obsid_lstfile):
			obsid = line.split()[0]

			# check gpg 
			flag_have_clevt = False			
			tmp_clevt_lst = glob.glob('%s/%s/xti/event_cl/ni%s_0mpu7_cl.evt*' % (self.param["input_data_directory"],obsid,obsid))
			for file in tmp_clevt_lst:
				if os.path.splitext(file)[-1] in ['.gz','.evt']:
					flag_have_clevt = True
					print('...find obsid %s cleaned event file.' % obsid)
					break 

			if obsid in ['#','%']:
				print('...comment out and skip the line: %s' % line)
				continue 
			elif obsid == 'exit' or obsid == 'quit' or obsid == 'exit()' or obsid == 'quit()' :
				print('...exit/quit of the obsid list.')
				break
			elif len(glob.glob('%s/%s' % (self.param["input_data_directory"],obsid))) == 0:
				print('...no directory and skip obsid ("input_data_directory" parameter): %s' % obsid)
				continue					
			elif not flag_have_clevt:
				print('...no clenaed event (gpg?) and skip obsid: %s' % obsid)
				continue					
			else:
				indir = glob.glob('%s/%s' % (self.param["input_data_directory"],obsid))[0]
				outdir = '%s/proc/%s' % (self.param['output_directory'],obsid)
				print('...set obsid:%s' % indir)
				self.nicerobs_lst.append(
					NicerObsID(obsid=obsid,indir=indir,outdir=outdir,param=self.param)
					)

		# --- initialization proclog ---
		self.proclog = NicerProcessLog('%s/nipipeline_proclog.csv' % self.param['output_directory'])
		for niobsid in self.nicerobs_lst:
			if str(niobsid.obsid) in list(self.proclog.df.index.tolist()):
				print('...ObsID %s is already in the existing list.' % niobsid.obsid)
				continue
			else:
				print('...ObsID %s is newly added to the existing list.' % niobsid.obsid)
				#self.proclog.df = self.proclog.df.append({'ObsID':str(niobsid.obsid)}, ignore_index=True)
				row = pd.Series(self.proclog.dict_null,name=str(niobsid.obsid))
				self.proclog.df = self.proclog.df.append(row)
				print(self.proclog.df)
		self.proclog.show()			

		if not os.path.exists(self.param["output_directory"]):
			cmd = 'mkdir -p %s' % self.param["output_directory"]
			print(cmd);os.system(cmd)
		self.proclog.show()			

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()

		print("*[NicerManager] %s DONE" % (sys._getframe().f_code.co_name))

	def make_directory(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_make_directory']:
			print("skip ... since flag_make_directory is %s" % self.param['flag_make_directory']) 
			return 0		

		self.proclog.df['mkdir'] = self.proclog.df['mkdir'].astype(str)
		for niobsid in self.nicerobs_lst:			
			flag = str(self.proclog.df.at[niobsid.obsid,'mkdir'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			try:
				niobsid.make_directory()
				self.proclog.df.at[niobsid.obsid,'mkdir'] = 'Done'		
				self.proclog.df.at[niobsid.obsid,'dir'] = '<a href="./%s">dir</a>' % (niobsid.outdir.replace(self.param['output_directory'],''))
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'dir'] = 'Error'
		
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()
		print("*[NicerManager] %s DONE" % (sys._getframe().f_code.co_name))

	def run_nicerl2(self):
		print("*[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nicerl2']:
			print("skip ... since flag_run_nicerl2 is %s" % self.param['flag_run_nicerl2']) 
			return 0
		
		self.proclog.df['nicerl2'] = self.proclog.df['nicerl2'].astype(str)
		self.proclog.df['object'] = self.proclog.df['object'].astype(str)
		self.proclog.df['date_obs'] = self.proclog.df['date_obs'].astype(str)
		self.proclog.df['exp(ks)'] = self.proclog.df['exp(ks)'].astype(str)
		self.proclog.df['#evt'] = self.proclog.df['#evt'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'nicerl2'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			try:
				niobsid.run_nicerl2()
				self.proclog.df.at[niobsid.obsid,'nicerl2'] = 'Done'

				self.proclog.df.at[niobsid.obsid,'object'] = niobsid.object
				self.proclog.df.at[niobsid.obsid,'date_obs'] = niobsid.date_obs
				self.proclog.df.at[niobsid.obsid,'exp(ks)'] = '%.2f' % (float(niobsid.exposure)*1e-3)
				self.proclog.df.at[niobsid.obsid,'#evt'] = '%.1e' % (float(niobsid.num_of_clevt))				
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'nicerl2'] = 'Error'				
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()
		print("*[NicerManager] %s DONE" % (sys._getframe().f_code.co_name))

	def run_nibackgen3C50(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nibackgen3C50']:
			print("skip ... since flag_run_nibackgen3C50 is %s" % self.param['flag_run_nibackgen3C50']) 
			return 0	

		self.proclog.df['spec'] = self.proclog.df['spec'].astype(str)
		self.proclog.df['3c50'] = self.proclog.df['3c50'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'3c50'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'3c50'] = 'Skip'	
				self.proclog.df.at[niobsid.obsid,'spec'] = 'Skip'
				continue 

			try:
				niobsid.run_nibackgen3C50()
				self.proclog.df.at[niobsid.obsid,'3c50'] = 'Done'	
				self.proclog.df.at[niobsid.obsid,'spec'] = '<a href="./%s">pdf</a>' % (niobsid.spec_pdf.replace(self.param['output_directory'],''))
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'3c50'] = 'Error'		
				self.proclog.df.at[niobsid.obsid,'spec'] = 'Error'
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()
		print("*[NicerManager] %s DONE" % (sys._getframe().f_code.co_name))

	def fit_avespec(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_fit_avespec']:
			print("skip ... since flag_run_nibackgen3C50 is %s" % self.param['flag_fit_avespec']) 
			return 0	

		self.proclog.df['xspec'] = self.proclog.df['xspec'].astype(str)
		self.proclog.df['fit'] = self.proclog.df['fit'].astype(str)
		self.proclog.df['rate'] = self.proclog.df['rate'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'fit'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'xspec'] = 'Skip'	
				self.proclog.df.at[niobsid.obsid,'fit'] = 'Skip'
				self.proclog.df.at[niobsid.obsid,'rate'] = 'Skip'
				continue 

			try:
				niobsid.fit_avespec()
				self.proclog.df.at[niobsid.obsid,'xspec'] = 'Done'	
				self.proclog.df.at[niobsid.obsid,'fit'] = '<a href="./%s">pdf</a>' % (niobsid.fit_avespec_pdf.replace(self.param['output_directory'],''))
				if float(niobsid.rate) < 1.0:				
					self.proclog.df.at[niobsid.obsid,'rate'] = '%.3f+/-%.3f' % (niobsid.rate,niobsid.rate_error)
				elif float(niobsid.rate) > 1000.0:	
					self.proclog.df.at[niobsid.obsid,'rate'] = '%.1e+/-%.1e' % (niobsid.rate,niobsid.rate_error)
				else:					
					self.proclog.df.at[niobsid.obsid,'rate'] = '%.1f+/-%.1f' % (niobsid.rate,niobsid.rate_error)					
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'xspec'] = 'Error'		
				self.proclog.df.at[niobsid.obsid,'fit'] = 'Error'
				self.proclog.df['rate'] = 'Error'
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()

	def plot_lightcurve(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_plot_lightcurve']:
			print("skip ... since flag_plot_lightcurve is %s" % self.param['flag_plot_lightcurve']) 
			return 0	

		self.proclog.df['lcurve'] = self.proclog.df['lcurve'].astype(str)
		self.proclog.df['lc'] = self.proclog.df['lc'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'lcurve'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'lcurve'] = 'Skip'	
				self.proclog.df.at[niobsid.obsid,'lc'] = 'Skip'
				continue 

			try:
				niobsid.plot_lightcurve()
				self.proclog.df.at[niobsid.obsid,'lcurve'] = 'Done'
				self.proclog.df.at[niobsid.obsid,'lc'] = '<a href="./%s">pdf</a>' % (niobsid.lc_pdf.replace(self.param['output_directory'],''))				
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'lcurve'] = 'Error'	
				self.proclog.df.at[niobsid.obsid,'lc'] = 'Error'
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()			

	def run_barycorr(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_barycorr']:
			print("skip ... since flag_devide_to_segment is %s" % self.param['flag_run_barycorr']) 
			return 0		

		self.proclog.df['bary'] = self.proclog.df['bary'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'bary'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'bary'] = 'Skip'	
				continue 

			try:
				niobsid.run_barycorr()
				self.proclog.df.at[niobsid.obsid,'bary'] = 'Done'
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'bary'] = 'Error'	

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()	

	def devide_to_segment(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_devide_to_segment']:
			print("skip ... since flag_devide_to_segment is %s" % self.param['flag_devide_to_segment']) 
			return 0		

		self.proclog.df['div2seg'] = self.proclog.df['div2seg'].astype(str)
		self.proclog.df['seg'] = self.proclog.df['seg'].astype(str)		
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'div2seg'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'div2seg'] = 'Skip'	
				self.proclog.df.at[niobsid.obsid,'seg'] = 'Skip'
				continue 

			try:
				niobsid.devide_to_segment()
				self.proclog.df.at[niobsid.obsid,'div2seg'] = 'Done'
				self.proclog.df.at[niobsid.obsid,'seg'] = '%d (<a href="./%s">dir</a>)' % (niobsid.number_of_segment, niobsid.dir_segment_main.replace(self.param['output_directory'],''))
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'div2seg'] = 'Error'	
				self.proclog.df.at[niobsid.obsid,'seg'] = 'Error'

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()	

	def fit_of_segment(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_fit_of_segment']:
			print("skip ... since fit_of_segment is %s" % self.param['flag_fit_of_segment']) 
			return 0		

		self.proclog.df['fit2seg'] = self.proclog.df['fit2seg'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'fit2seg'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'fit2seg'] = 'Skip'	
				continue 

			try:
				niobsid.fit_of_segment()
				self.proclog.df.at[niobsid.obsid,'fit2seg'] = 'Done'
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)				
				self.proclog.df.at[niobsid.obsid,'fit2seg'] = 'Error'	

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()	

	def devide_to_block(self):
		print("-[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_devide_to_block']:
			print("skip ... since flag_devide_to_block is %s" % self.param['flag_devide_to_block']) 
			return 0		

		self.proclog.df['div2block'] = self.proclog.df['div2block'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'div2block'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'div2block'] = 'Skip'	
				continue 

			try:
				niobsid.devide_to_block()
				self.proclog.df.at[niobsid.obsid,'div2block'] = 'Done'
				#self.proclog.df.at[niobsid.obsid,'seg'] = '%d (<a href="./%s">dir</a>)' % (niobsid.number_of_segment, niobsid.dir_segment_main.replace(self.param['output_directory'],''))
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)
				self.proclog.df.at[niobsid.obsid,'div2block'] = 'Error'	
				#self.proclog.df.at[niobsid.obsid,'seg'] = 'Error'

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()	


	def fit_of_block(self):
		print("-[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_fit_of_block']:
			print("skip ... since flag_devide_to_block is %s" % self.param['flag_fit_of_block']) 
			return 0		

		self.proclog.df['#block'] = self.proclog.df['#block'].astype(str)
		self.proclog.df['fit2block'] = self.proclog.df['fit2block'].astype(str)
		for niobsid in self.nicerobs_lst:
			flag = str(self.proclog.df.at[niobsid.obsid,'fit2block'])
			if self.param['flag_process_overwrite'] == False and flag == 'True':
				print("...skip, already processed.")
				continue

			if not niobsid.flag_clevt_has_events:
				self.proclog.df.at[niobsid.obsid,'fit2block'] = 'Skip'	
				self.proclog.df.at[niobsid.obsid,'#block'] =  'Skip'

				continue 

			try:
				niobsid.fit_of_block()
				self.proclog.df.at[niobsid.obsid,'fit2block'] = 'Done'
				self.proclog.df.at[niobsid.obsid,'#block'] =  '%d' % (niobsid.total_number_of_block)
			except:
				print("Exception: %s" % sys._getframe().f_code.co_name)
				self.proclog.df.at[niobsid.obsid,'fit2block'] = 'Error'	
				self.proclog.df.at[niobsid.obsid,'#block'] =  'Error'

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()	

