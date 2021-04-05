# -*- coding: utf-8 -*-

import os 
import sys
import glob
import yaml 
import pandas as pd
import inspect
from time import sleep

from astropy.io import fits
from astropy.time import Time

def nibackgen3C50(outdir,basename,totspec_name,bkgspec_name,
		bkgidxdir,bkglibdir,gainepoch,
		mode='obsid',obsid=None,clfile=None,ufafile=None,
		dtmin=20.0,dtmax=60.6,hbgcut=0.5,rmf=None,arf=None):
	print("\n[func:%s]" % sys._getframe().f_code.co_name)

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

	print("\n[func:%s] DONE" % sys._getframe().f_code.co_name)
	return 0

def plot_xspec(totspec,bkgspec,
	xspec_rebin_sigma=5,xspec_rebin_maxnum=100,
	xspec_emin=0.2,xspec_emax=10.0,xspec_plot_type='ld',
	xspec_ymin=1e-4,xspec_ymax=100):
	print("\n[func:%s]" % sys._getframe().f_code.co_name)

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

	print("\n[func:%s] Done" % sys._getframe().f_code.co_name)

def plot_lcurve(eventfile,outdir,lc_energy_bands,lc_time_bin_sec):
	print("\n[func:%s]" % sys._getframe().f_code.co_name)

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
	dump += 'plot=yes plotdev="/xw" <<EOF\n'
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

	print("\n[func:%s] Done" % sys._getframe().f_code.co_name)

class NicerObsID():
	def __init__(self,obsid,indir,outdir,param):
		self.param = param; 
		self.obsid = obsid; 
		self.indir = indir; 
		self.outdir = outdir; 
		print("[NicerObsID %s] A new NicerObs is generated" % self.obsid)

		#self.errorlog = '%s/errorlog.txt' % self.param['output_directory']
		self.setup_yamlfile = '%s/ni%s_result.yaml' % (self.outdir,self.obsid)
		self.have_event_files()
		sleep(1) # this is needed, otherwise sometimes flag_clevt_has_events can not be ste.

	def have_event_files(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		self.clevt = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (self.outdir,self.obsid)		
		if os.path.exists(self.clevt):
			self.flag_exists_clevt = True
			hdu = fits.open(self.clevt)
			if len(hdu['EVENTS'].data) > 0:
				self.flag_clevt_has_events = True
			else:
				self.flag_clevt_has_events = False
		else:
			self.flag_exists_clevt = False

		self.ufaevt = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt' % (self.outdir,self.obsid)		
		if os.path.exists(self.ufaevt):
			self.flag_exists_ufaevt = True
		else:
			self.flag_exists_ufaevt = False
		print("\n[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0 

	def make_directory(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		cmd  = 'rm -rf %s;\n' % (self.outdir)
		cmd += 'mkdir -p %s/xti;\n' % (self.outdir)
		cmd += 'ln -s %s/{auxil,log} %s/;\n' % (self.indir,self.outdir)
		cmd += 'ln -s %s/xti/{event_uf,hk} %s/xti/;\n' % (self.indir,self.outdir)	
		print(cmd);os.system(cmd)
		self.dump_setup_to_yamlfile()
		print("\n[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0 

	def dump_setup_to_yamlfile(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

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
		print("\n[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def run_nicerl2(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
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
		self.dump_setup_to_yamlfile()	
		print("\n[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def run_nibackgen3C50(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		fversion = os.popen('fversion').read().rstrip()		

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

		totspec = '%s/spec/ni%s_3c50_ave_tot.pi' % (self.outdir,self.obsid)
		bkgspec = '%s/spec/ni%s_3c50_ave_bkg.pi' % (self.outdir,self.obsid)		
		plot_xspec(
			totspec=totspec,
			bkgspec=bkgspec,
			xspec_rebin_sigma=self.param['xspec_rebin_sigma'],
			xspec_rebin_maxnum=self.param['xspec_rebin_maxnum'],
			xspec_emin=self.param['xspec_emin'],
			xspec_emax=self.param['xspec_emax'],
			xspec_plot_type=self.param['xspec_plot_type'],
			xspec_ymin=self.param['xspec_ymin'],
			xspec_ymax=self.param['xspec_ymax'])
		self.spec_pdf = '%s/spec/ni%s_3c50_ave_tot.pdf' % (self.outdir,self.obsid)
		self.dump_setup_to_yamlfile()
		print("\n[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0

	def plot_lightcurve(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		plot_lcurve(
			eventfile=self.clevt,
			outdir=self.outdir,
			lc_energy_bands=self.param['lc_energy_bands'],
			lc_time_bin_sec=self.param['lc_time_bin_sec']
			)

		self.lc_pdf = '%s/lc/ni%s_0mpu7_cl_ene.pdf' % (self.outdir,self.obsid)

		self.dump_setup_to_yamlfile()

		print("\n[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0		

	def devide_to_segment(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		suboutdir = '%s/segment' % self.outdir
		cmd  = 'rm -rf %s;\n' % suboutdir
		cmd += 'mkdir -p %s\n' % suboutdir
		print(cmd);os.system(cmd)

		# re-arrange GTIs
		cmd = 'rm -f tmp_gti.txt; group_gti.py %s tmp_gti.txt --gti_separation_threshold %.2f' % (self.clevt,self.param['gti_separation_threshold'])
		print(cmd);os.system(cmd)

		# convert text file GTI format to fits file.
		segment_gti = '%s/ni%s_segment.gti' % (suboutdir,self.obsid)
		cmd = 'fconv_txt2gti.py -i tmp_gti.txt -o %s; rm -f tmp_gti.txt' % segment_gti 
		print(cmd);os.system(cmd)

		hdu = fits.open(segment_gti)
		number_of_segment = len(hdu['GTI'].data)

		self.nigsegment_list = []
		for i in range(number_of_segment):
			segment_num = i + 1 
			segment_dir = '%s/segment%003d' % (suboutdir,segment_num)
			print("----- segment_dir = %s" % segment_dir)
			cmd  = 'mkdir -p %s;\n' % segment_dir
			print(cmd);os.system(cmd)

			hdu = fits.open(segment_gti)
			hdu['GTI'].data = hdu['GTI'].data[[i]]
			hdu.writeto('%s/tmp_sel.gti' % segment_dir)

			segment_basename = 'ni%s_0mpu7_cl_seg%03d' % (self.obsid,segment_num)
			segment_fitsfile = '%s.evt' % segment_basename
			cmd  = 'xselect_gtifilter.py '
			cmd += '-i %s ' % self.clevt
			cmd += '-g %s/tmp_sel.gti ' % segment_dir
			cmd += '-o %s/%s;\n' % (segment_dir,segment_fitsfile)
			print(cmd);os.system(cmd)

			segment_ufa_basename = 'ni%s_0mpu7_ufa_seg%03d' % (self.obsid,segment_num)
			segment_ufa_fitsfile = '%s.evt' % segment_ufa_basename
			cmd  = 'xselect_gtifilter.py '
			cmd += '-i %s ' % self.ufaevt
			cmd += '-g %s/tmp_sel.gti ' % segment_dir
			cmd += '-o %s/%s;\n' % (segment_dir,segment_ufa_fitsfile)
			print(cmd);os.system(cmd)

			nisegment = NicerSegment(self,
				segment_num,segment_dir,segment_basename)
			if len(nisegment.segment_clhdu['EVENTS'].data) == 0:
				print("... No event. skip.")
				continue			

			if nisegment.nbint > 0: 
				nisegment.plot_lightcurve()
			nisegment.run_nibackgen3C50()
			self.nigsegment_list.append(nisegment)

		"""
		self.df_summary = pd.DataFrame(columns=['ObsID','segment','lc','spec'])
		#for i in range(number_of_segment):
		for nisegment in self.nigsegment_list:
			if nisegment.lc_pdf == None:
				lc_link = 'N/A'
			else:
				lc_link = '<a href="./%s">pdf</a>' % nisegment.lc_pdf.replace(suboutdir,'')
			if nisegment.spec_pdf == None:
				spec_link = 'N/A'
			else:
				spec_link = '<a href="./%s">pdf</a>' % nisegment.spec_pdf.replace(suboutdir,'')			
			self.df_summary = self.df_summary.append({'ObsID':self.obsid,
				'segment':nisegment.segment_num,
				'lc':lc_link,'spec':spec_link}, ignore_index=True)

		self.show_dataframe_summay()
		self.write_dataframe_summary()
		"""
		print("\n[NicerObsID %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0		

class NicerSegment():
	def __init__(self,parent_nicerobsid,
			segment_num,segment_dir,segment_basename):
		self.parent_nicerobsid = parent_nicerobsid
		self.segment_num = segment_num; 
		self.segment_dir = segment_dir
		self.segment_basename = segment_basename
		print("\n[NicerSegment %d] A new NicerSegment is generated from a parent NicerObdID %s" % (
			self.segment_num,self.parent_nicerobsid.obsid))

		self.obsid = self.parent_nicerobsid.obsid
		self.param = self.parent_nicerobsid.param

		self.segment_gti = '%s/%s.gti' % (self.segment_dir,self.segment_basename)
		self.segment_clevt = '%s/%s.evt' % (self.segment_dir,self.segment_basename)
		self.segment_ufaevt = '%s/ni%s_0mpu7_ufa_seg%03d.evt' % (self.segment_dir,self.obsid,self.segment_num)

		self.segment_clhdu = fits.open(self.segment_clevt)

		self.totspec_name = 'ni%s_seg%03d_3c50_tot' % (self.obsid,self.segment_num)
		self.bkgspec_name = 'ni%s_seg%03d_3c50_bkg' % (self.obsid,self.segment_num)
		self.totspec = '%s/%s.pi' % (self.segment_dir,self.totspec_name)
		self.bkgspec = '%s/%s.pi' % (self.segment_dir,self.bkgspec_name)

		self.spec_pdf = None
		self.spec_link = None
		self.lc_pdf = None
		self.lc_link = None		

		self.set_title()		

	def set_title(self):
		print("\n[NicerSegment %s/%s] %s" % (self.segment_num,self.obsid,sys._getframe().f_code.co_name))

		self.tstart = float(self.segment_clhdu['GTI'].data[0][0])
		self.tstop = float(self.segment_clhdu['GTI'].data[-1][1])
		self.duration = self.tstop - self.tstart

		hdu = fits.open(self.segment_clevt)
		self.exposure = float(hdu[0].header['EXPOSURE'])

		self.nbint = round(self.duration/float(self.param['lc_time_bin_sec']))
		self.title = 'ObsID:%s GTI:%003d MJD:%.2f-%.2f (%.1f s)' % (
			self.obsid,self.segment_num,self.tstart,self.tstop,self.exposure)
		print(self.duration,self.nbint,self.param['lc_time_bin_sec'])
		print("\n[NicerSegment %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))

	def plot_lightcurve(self):
		print("\n[NicerSegment %s/%s] %s" % (self.segment_num,self.obsid,sys._getframe().f_code.co_name))

		plot_lcurve(
			eventfile=self.segment_clevt,
			outdir=self.segment_dir,
			lc_energy_bands=self.param['lc_energy_bands'],
			lc_time_bin_sec=self.param['lc_time_bin_sec']
			)

		self.lc_pdf = '%s/%s_ene.pdf' % (self.segment_dir,self.segment_basename)

		print("\n[NicerSegment %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0		

	def run_nibackgen3C50(self):
		print("\n[NicerSegment %s/%s] %s" % (self.segment_num,self.obsid,sys._getframe().f_code.co_name))

		flag_nibackgen3c50 = nibackgen3C50(
			outdir=self.segment_dir,
			basename=self.segment_basename,
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
			clfile=self.segment_clevt,
			ufafile=self.segment_ufaevt,
			)

		if flag_nibackgen3c50 == 0:
			totspec = '%s/spec/%s.pi' % (self.segment_dir,self.totspec_name)
			bkgspec = '%s/spec/%s.pi' % (self.segment_dir,self.bkgspec_name)		
			plot_xspec(
				totspec=totspec,
				bkgspec=bkgspec,
				xspec_rebin_sigma=self.param['xspec_rebin_sigma'],
				xspec_rebin_maxnum=self.param['xspec_rebin_maxnum'],
				xspec_emin=self.param['xspec_emin'],
				xspec_emax=self.param['xspec_emax'],
				xspec_plot_type=self.param['xspec_plot_type'],
				xspec_ymin=self.param['xspec_ymin'],
				xspec_ymax=self.param['xspec_ymax'])
			self.spec_pdf = '%s/spec/%s.pdf' % (self.segment_dir,self.totspec_name)
		#self.dump_setup_to_yamlfile()

		print("\n[NicerSegment %s] %s DONE" % (self.obsid,sys._getframe().f_code.co_name))
		return 0	


class NicerProcessLog():
	def __init__(self,csvfile):
		print("\n[NicerProcessLog] A new NicerProcessLog object is generated.")

		self.csvfile = csvfile
		self.htmlfile = os.path.splitext(self.csvfile)[0] + '.html'

		if os.path.exists(self.csvfile):
			print("...reading %s." % self.csvfile )
			self.df = pd.read_csv(csvfile,index_col=0)
			#self.df.dtype=str
			self.df = self.df.astype(str)
		else:
			print("...making a new process log.")
			self.df = pd.DataFrame(
				columns=['ObsID','dir','nicerl2','3c50','lcurve','segment','directory','spec','lc'],
				dtype=str,
				)
		self.show()

	def show(self):
		print("******* summary ******")
		print(self.df)
		print("**********************")		

	def write_to_csvfile(self):
		print("\n[NicerProcessLog] %s" % (sys._getframe().f_code.co_name))
		self.df.to_csv(self.csvfile)		

	def write_to_htmlfile(self):
		print("\n[NicerProcessLog] %s" % (sys._getframe().f_code.co_name))
		raw_html = self.df.to_html(open(self.htmlfile,'w'))
		cmd = 'sed -e "s/&lt;/</g" %s > tmp.html' % self.htmlfile
		print(cmd);os.system(cmd)
		cmd = 'sed -e "s/&gt;/>/g" tmp.html > tmp2.html' 
		print(cmd);os.system(cmd)		
		cmd = 'mv tmp2.html %s; rm -f tmp.html tmp2.html' % self.htmlfile
		print(cmd);os.system(cmd)		


class NicerManager():
	def __init__(self,setup_yamlfile,obsid_lstfile):
		print("\n[NicerManager] A new NicerManager object is generated.")

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
		print(self.proclog.df['ObsID'])
		for niobsid in self.nicerobs_lst:
			if str(niobsid.obsid) in list(self.proclog.df['ObsID']):
				print('...ObsID %s is already in the existing list.' % niobsid.obsid)
				continue
			else:
				print('...ObsID %s is newly added to the existing list.' % niobsid.obsid)
				self.proclog.df = self.proclog.df.append({'ObsID':str(niobsid.obsid)}, ignore_index=True)
		self.proclog.show()			

		if not os.path.exists(self.param["output_directory"]):
			cmd = 'mkdir -p %s' % self.param["output_directory"]
			print(cmd);os.system(cmd)
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()

	def make_directory(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_make_directory']:
			print("skip ... since flag_make_directory is %s" % self.param['flag_make_directory']) 
			return 0		
		dict_flag = {}
		dict_output = {}
		for niobsid in self.nicerobs_lst:
			flag_dir = str(self.proclog.df[self.proclog.df['ObsID'] == niobsid.obsid]['dir'][0])
			if self.param['flag_process_overwrite'] == False and flag_dir == 'True':
				print("...skip, already processed.")
				continue
			try:
				niobsid.make_directory()
				flag = True
				output = '<a href="./%s">dir</a>' % (niobsid.outdir.replace(self.param['output_directory'],''))
			except:
				flag = False
				output = 'ERROR'
			dict_flag[flag] = niobsid.obsid
			dict_output[output] = niobsid.obsid
		if len(dict_flag) > 0:
			self.proclog.df['dir'] = dict_flag			
			self.proclog.df['directory'] = dict_output
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()

	def run_nicerl2(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nicerl2']:
			print("skip ... since flag_run_nicerl2 is %s" % self.param['flag_run_nicerl2']) 
			return 0
		dict_flag = {}						
		for niobsid in self.nicerobs_lst:
			flag_nicerl2 = str(self.proclog.df[self.proclog.df['ObsID'] == niobsid.obsid]['nicerl2'][0])
			if self.param['flag_process_overwrite'] == False and flag_nicerl2 == 'True':
				print("...skip, already processed.")
				continue

			try:
				niobsid.run_nicerl2()
				flag = True
			except:
				flag = False
			dict_flag[flag] = niobsid.obsid
		if len(dict_flag) > 0:
			self.proclog.df['nicerl2'] = dict_flag
		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()

	def run_nibackgen3C50(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nibackgen3C50']:
			print("skip ... since flag_run_nibackgen3C50 is %s" % self.param['flag_run_nibackgen3C50']) 
			return 0	

		dict_flag = {}		
		dict_output = {}		
		for niobsid in self.nicerobs_lst:
			flag_spec = str(self.proclog.df[self.proclog.df['ObsID'] == niobsid.obsid]['3c50'][0])
			if self.param['flag_process_overwrite'] == False and flag_spec == 'True':
				print("...skip, already processed.")
				continue

			try:
				niobsid.run_nibackgen3C50()
				flag = True
				output = '<a href="./%s">dir</a>' % (niobsid.spec_pdf.replace(self.param['output_directory'],''))				
			except:
				flag = False
				output = 'ERROR'
			dict_flag[flag] = niobsid.obsid
			dict_output[output] = niobsid.obsid

		if len(dict_flag) > 0:
			self.proclog.df['3c50'] = dict_flag
			self.proclog.df['spec'] = dict_output

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()

	def plot_lightcurve(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_plot_lightcurve']:
			print("skip ... since flag_plot_lightcurve is %s" % self.param['flag_plot_lightcurve']) 
			return 0	

		dict_flag = {}		
		dict_output = {}		
		for niobsid in self.nicerobs_lst:
			flag_lc = str(self.proclog.df[self.proclog.df['ObsID'] == niobsid.obsid]['lc'][0])
			if self.param['flag_process_overwrite'] == False and flag_lc == 'True':
				print("...skip, already processed.")
				continue

			try:
				niobsid.plot_lightcurve()
				flag = True
				output = '<a href="./%s">dir</a>' % (niobsid.lc_pdf.replace(self.param['output_directory'],''))				
			except:
				flag = False
				output = 'ERROR'
			dict_flag[flag] = niobsid.obsid
			dict_output[output] = niobsid.obsid

		if len(dict_flag) > 0:
			self.proclog.df['lcurve'] = dict_flag
			self.proclog.df['lc'] = dict_output

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()			

	def devide_to_segment(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_devide_to_segment']:
			print("skip ... since flag_devide_to_segment is %s" % self.param['flag_devide_to_segment']) 
			return 0		

		dict_flag = {}		
		dict_output = {}		
		for niobsid in self.nicerobs_lst:
			flag_segment = str(self.proclog.df[self.proclog.df['ObsID'] == niobsid.obsid]['segment'][0])
			if self.param['flag_process_overwrite'] == False and flag_segment == 'True':
				print("...skip, already processed.")
				continue

			niobsid.devide_to_segment()
		"""
			try:
				niobsid.devide_to_segment()
				flag = True
				#output = '<a href="./%s">dir</a>' % (niobsid.lc_pdf.replace(self.param['output_directory'],''))				
			except:
				flag = False
				#output = 'ERROR'
			dict_flag[flag] = niobsid.obsid
			#dict_output[output] = niobsid.obsid

		if len(dict_flag) > 0:
			self.proclog.df['segment'] = dict_flag
			#self.proclog.df['lc'] = dict_output

		self.proclog.write_to_csvfile()
		self.proclog.write_to_htmlfile()	
		"""


"""
class NicerObsID
	def run_barycorr(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		suboutdir = '%s/bary' % self.outdir
		cmd  = 'rm -rf %s;\n' % suboutdir
		cmd += 'mkdir -p %s\n' % suboutdir
		print(cmd);os.system(cmd)

		if not self.flag_clevt_has_events:
			flog = '%s/barycorr_%s.log' % (suboutdir,self.obsid)
			f = open(flog,'w')
			message = 'no cleaned event in the original file %s' % self.clevt
			print(message)
			f.write(message)
			f.close()
			return -1 

		clevt = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (self.outdir,self.obsid)
		hdu = fits.open(clevt)
		if self.param['barycorr_RA_deg'] == "None":
			self.param['barycorr_RA_deg'] = float(hdu['EVENTS'].header['RA_OBJ'])
		if self.param['barycorr_DEC_deg'] == "None":
			self.param['barycorr_DEC_deg'] = float(hdu['EVENTS'].header['DEC_OBJ'])
		print(self.param['barycorr_RA_deg'])
		print(self.param['barycorr_DEC_deg'])		

		baryevt = '%s/ni%s_0mpu7_cl_bary.evt' % (suboutdir,self.obsid)
		fcmd = '%s/barycorr_%s.sh' % (suboutdir,self.obsid)
		flog = '%s/barycorr_%s.log' % (suboutdir,self.obsid)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'
		dump += 'barycorr infile=%s ' % clevt
		dump += 'outfile=%s ' % baryevt
		dump += 'orbitfiles="%s/auxil/ni%s.orb.gz" ' % (self.outdir,self.obsid)
		dump += 'ra=%.6f dec=%.6f ' % (self.param['barycorr_RA_deg'],self.param['barycorr_DEC_deg'])
		dump += 'refframe=ICRS ephem=JPLEPH.430 '
		dump += '> %s 2>&1 ' % flog
		dump += '\n'
		f.write(dump)
		f.close()	

		# run the script.
		cmd  = 'chmod +x %s\n' % fcmd
		cmd += './%s' % fcmd
		print(cmd);os.system(cmd)	

		self.dump_setup_to_yamlfile()


	def fit_of_segment(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		suboutdir = '%s/segment' % self.outdir
		segment_gti = '%s/ni%s_segment.gti' % (suboutdir,self.obsid)
		hdu = fits.open(segment_gti)
		number_of_segment = len(hdu['GTI'].data)

		for i in range(number_of_segment):
			segment_num = i + 1 
			segment_dir = '%s/ep%003d' % (suboutdir,segment_num)
			print("----- segment_dir = %s" % segment_dir)
			segment_basename = 'ni%s_0mpu7_cl_ep%03d' % (self.obsid,segment_num)
			segment_fitsfile = '%s.evt' % segment_basename
			nisegment = NicerSegment(self,segment_num,segment_dir,segment_basename)
			if len(nisegment.segment_clhdu['EVENTS'].data) == 0:
				print("... No event. skip.")
				continue
			nisegment.set_title()
			try:
				nisegment.fit_spectrum()
			except:				
				message  = "Error: [NicerObsID %s] %s\n" % (self.obsid,sys._getframe().f_code.co_name)
				message += "        self.obsid=%s, segment_num=%d\n" % (self.obsid, segment_num)
				f = open(self.errorlog,'w')
				f.write(message)
				f.close()

		dict_append = {}				
		for nisegment in self.nigsegment_list:
			link = '<a href="./%s">pdf</a>' % nisegment.fit_pdf.replace(suboutdir,'')
			dict_append[link] = nisegment.segment_num
		self.df_summary['fit'] = dict_append

		self.show_dataframe_summay()
		self.write_dataframe_summary()			

	def show_dataframe_summay(self):
		print("******* summary ******")
		print(self.df_summary)
		print("**********************")		

	def write_dataframe_summary(self):
		self.summary_csv = '%s/segment/ni%s_segment.csv' % (self.outdir,self.obsid)
		self.summary_html = '%s/segment/ni%s_segment.html' % (self.outdir,self.obsid)
		self.df_summary.to_csv(self.summary_csv)

		raw_html = self.df_summary.to_html(open(self.summary_html,'w'))
		cmd = 'sed -e "s/&lt;/</g" %s > tmp.html' % self.summary_html
		print(cmd);os.system(cmd)
		cmd = 'sed -e "s/&gt;/>/g" tmp.html > tmp2.html' 
		print(cmd);os.system(cmd)		
		cmd = 'mv tmp2.html %s; rm -f tmp.html tmp2.html' % self.summary_html
		print(cmd);os.system(cmd)		



clss segment 


	def fit_spectrum(self):
		cmd  = 'xspec_fit.py '
		cmd += '%s ' % self.totspec
		cmd += '-o %s/fit ' % self.segment_dir
		cmd += '-b %s ' % self.bkgspec
		cmd += '-r %s ' % self.param['xspec_rmf']
		cmd += '-a %s ' % self.param['xspec_arf']
		cmd += '-m %s ' % self.param['xspec_model']
		cmd += '-s %d -n %d ' % (self.param['xspec_rebin_sigma'],self.param['xspec_rebin_maxnum'])
		cmd += '--fitemin %.1f --fitemax %.1f ' % (self.param['xspec_emin'],self.param['xspec_emax'])
		cmd += '--rateband %s ' % self.param['xspec_rateband']
		cmd += '--fluxband %s \n' % self.param['xspec_fluxband']
		print(cmd);os.system(cmd)

		self.fit_pdf = '%s/fit/%s.pdf' % (self.segment_dir,fspec)





	def run_barycorr(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_barycorr']:
			print("skip ... since flag_run_barycorr is %s" % self.param['flag_run_barycorr']) 
			return 0			
		for niobsid in self.nicerobs_lst:
			niobsid.run_barycorr()

	def fit_of_segment(self):
		print("\n[NicerManager] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_fit_of_segment']:
			print("skip ... since flag_fit_of_segment is %s" % self.param['flag_fit_of_segment']) 
			return 0		
		for niobsid in self.nicerobs_lst:
			niobsid.fit_of_segment()
"""			
