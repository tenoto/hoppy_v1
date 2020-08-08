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

class NicerObsID():
	def __init__(self,obsid,indir,outdir,param):
		self.param = param; 
		self.obsid = obsid; 
		self.indir = indir; 
		self.outdir = outdir; 
		print("[NicerObsID %s] A new NicerObs is generated" % self.obsid)

		self.setup_yamlfile = '%s/ni%s_result.yaml' % (self.outdir,self.obsid)
		self.have_event_files()
		sleep(1) # this is needed, otherwise sometimes flag_clevt_has_events can not be ste.

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

	def have_event_files(self):
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

	def make_directory(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		cmd  = 'rm -rf %s;\n' % (self.outdir)
		cmd += 'mkdir -p %s/xti;\n' % (self.outdir)
		cmd += 'ln -s %s/{auxil,log} %s/;\n' % (self.indir,self.outdir)
		cmd += 'ln -s %s/xti/{event_uf,hk} %s/xti/;\n' % (self.indir,self.outdir)	
		print(cmd);os.system(cmd)
		self.dump_setup_to_yamlfile()

	def run_nicerl2(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		fversion = os.popen('fversion').read().rstrip()		
	
		# prepare a script for nicerl2 for each ObsID
		fcmd = '%s/nicerl2_%s.sh' % (self.outdir,self.obsid)
		flog = '%s/nicerl2_%s.log' % (self.outdir,self.obsid)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'
		dump += 'nicerl2 indir=%s ' % (self.outdir)
		dump += 'picalfile=%s ' % self.param['nicerl2_gcalfile']
		dump += 'clobber=yes '
		dump += '> %s 2>&1 ' % flog
		dump += '\n'
		f.write(dump)
		f.close()

		# run the script.
		cmd  = 'chmod +x %s\n' % fcmd
		cmd += './%s' % fcmd
		print(cmd)
		os.system(cmd)

		f = open(flog,'a')
		f.write(fversion+'\n')
		f.close()

		self.have_event_files()
		self.dump_setup_to_yamlfile()		

	def run_nibackgen3C50(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))
		fversion = os.popen('fversion').read().rstrip()		

		suboutdir = '%s/spec' % self.outdir
		cmd  = 'rm -rf %s;\n' % (suboutdir)
		cmd += 'mkdir -p %s\n' % (suboutdir)
		print(cmd);os.system(cmd)

		if not self.flag_clevt_has_events:
			flog = '%s/nibackgen3C50_%s.log' % (suboutdir,self.obsid)
			f = open(flog,'w')
			message = 'no cleaned event in the original file %s' % self.clevt
			print(message)
			f.write(message)
			f.close()
			return -1 

		totspec_name = 'ni%s_3c50_tot' % self.obsid
		bkgspec_name = 'ni%s_3c50_bkg' % self.obsid	

		# prepare a script for nicerl2 for each self.obsid
		fcmd = '%s/nibackgen3C50_%s.sh' % (suboutdir,self.obsid)
		flog = '%s/nibackgen3C50_%s.log' % (suboutdir,self.obsid)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'
		dump += 'nibackgen3C50 ' 
		dump += 'rootdir=\'%s\' ' % self.param['output_directory']
		dump += 'obsid=\'%s\' ' % self.obsid
		dump += 'bkgidxdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']
		dump += 'bkglibdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']	
		dump += 'gainepoch=\'%s\' ' % self.param['nibackgen3c50_gainepoch']
		dump += 'totspec=\'%s\' ' % (totspec_name)
		dump += 'bkgspec=\'%s\' ' % (bkgspec_name)
		dump += 'dtmin=%.1f dtmax=%.1f hbgcut=%.1f ' % (
			self.param['nibackgen3c50_dtmin'],
			self.param['nibackgen3c50_dtmax'],
			self.param['nibackgen3c50_hbgcut'])
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

		f = open(flog,'a')
		f.write(fversion+'\n')
		f.close()		

		if os.path.exists(totspec_name+'.pi'):
			self.flag_exists_totspec = True
		else:
			self.flag_exists_totspec = False
			print("spectral file %s was not generated." % totspec_name)
			return -1

		if os.path.exists(bkgspec_name+'.pi'):
			self.flag_exists_bkgspec = True
		else:
			self.flag_exists_bkgspec = False
			print("spectral file %s was not generated." % bkgspec_name)
			return -1			

		cmd = 'mv %s*.pi %s*.pi %s' % (totspec_name,bkgspec_name,suboutdir)
		print(cmd);os.system(cmd)		

		totspec = '%s/%s.pi' % (suboutdir,totspec_name)
		bkgspec = '%s/%s.pi' % (suboutdir,bkgspec_name)

		cmd  = 'fparkey %s %s RESPFILE\n' % (self.param['xspec_rmf'],totspec)
		cmd += 'fparkey %s %s ANCRFILE\n' % (self.param['xspec_arf'],totspec)
		#cmd += 'fparkey %s %s BACKFILE\n' % (bkgspec,totspec)
		cmd += 'fparkey %s %s RESPFILE\n' % (self.param['xspec_rmf'],bkgspec)
		cmd += 'fparkey %s %s ANCRFILE\n' % (self.param['xspec_arf'],bkgspec)
		print(cmd);os.system(cmd)

		hdu = fits.open(totspec)
		keyword_dateobs = hdu[0].header['DATE-OBS']
		keyword_object = hdu[0].header['OBJECT']		
		keyword_exposure = hdu[0].header['EXPOSURE']			
		title = '%s ObsID:%s (%s, %.1f sec)' % (keyword_object, self.obsid,keyword_dateobs,keyword_exposure)
		fspec = "ni%s_tot_spec" % self.obsid
		cmd  = "xspec <<EOF\n"
		cmd += "data 1 %s\n" % (totspec)
		cmd += "back 1 %s\n" % (bkgspec)
		cmd += "data 2 %s\n" % (bkgspec)
		cmd += "data 3 %s\n" % (totspec)	
		cmd += "back 3 none\n"
		cmd += "setplot energy\n"
		cmd += "setplot rebin %d %d\n" % (self.param['xspec_rebin_sigma'],self.param['xspec_rebin_maxnum'])
		cmd += "ignore 1-3:**-%.1f,%.1f-**\n" % (self.param['xspec_emin'],self.param['xspec_emax'])
		cmd += "iplot %s\n" % self.param['xspec_plot_type']
		cmd += "time off\n"
		cmd += "lwid 5\n"
		cmd += "lwid 5 on 1..3\n"
		cmd += "la t %s\n" % title
		cmd += "r x %.1f %.1f\n" % (self.param['xspec_emin'],self.param['xspec_emax'])
		#cmd += "r y 1e-4 100.0\n"
		cmd += "r y %.1f %.1f\n" % (self.param['xspec_ymin'],self.param['xspec_ymax'])
		cmd += "hard %s.ps/cps\n" % fspec
		cmd += "exit\n"
		cmd += "exit\n"
		cmd += "EOF\n"
		print(cmd)
		os.system(cmd)

		cmd  = "ps2pdf %s.ps\n" % fspec
		cmd += "rm -f %s.ps\n" % fspec
		cmd += "mv %s.pdf %s/\n" % (fspec,suboutdir)
		print(cmd)
		os.system(cmd)

		self.dump_setup_to_yamlfile()

	def plot_lightcurve(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		suboutdir = '%s/lc' % self.outdir
		cmd  = 'rm -rf %s;\n' % suboutdir
		cmd += 'mkdir -p %s\n' % suboutdir
		print(cmd);os.system(cmd)

		if not self.flag_clevt_has_events:
			flog = '%s/lcurve_%s.log' % (suboutdir,self.obsid)
			f = open(flog,'w')
			message = 'no cleaned event in the original file %s' % self.clevt
			print(message)
			f.write(message)
			f.close()
			return -1 

		for eband in self.param['lc_energy_bands']:
			emin,emax = self.param['lc_energy_bands'][eband]
			fenesel = 'ni%s_0mpu7_cl_%s.evt' % (self.obsid,eband)
			cmd += "fselect_energy.py %s/xti/event_cl/ni%s_0mpu7_cl.evt %s %.1f %.1f;\n" % (self.outdir,self.obsid,fenesel,emin,emax)
			cmd += "mv %s %s/lc;\n" % (fenesel,self.outdir)
		print(cmd)
		os.system(cmd)

		hdu = fits.open('%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (self.outdir,self.obsid))
		t = Time([hdu[1].header['DATE-OBS'],hdu[1].header['DATE-END']],format='isot', scale='utc')
		key_date_obs = hdu[1].header['DATE-OBS']      
		key_date_end = hdu[1].header['DATE-END']      
		key_exposure = hdu[1].header['EXPOSURE']
		key_target = hdu[1].header['OBJECT']
		key_mjd_start = t.mjd[0]
		key_mjd_stop = t.mjd[1]       
		key_onsource = (key_mjd_stop-key_mjd_start)*24.0*60.0*60.0
		nbint = round(float(key_onsource) / float(self.param['lc_time_bin_sec']))
		title = 'ObsID:%s %s MJD:%.2f (exposure %.1f s)' % (self.obsid,key_date_obs,key_mjd_start,key_exposure)

		fcmd = '%s/ni%s_0mpu7_cl_ene.sh' % (suboutdir,self.obsid)
		flog = '%s/ni%s_0mpu7_cl_ene.log' % (suboutdir,self.obsid)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'
		flcfile = '%s/ni%s_0mpu7_cl_ene.flc' % (suboutdir,self.obsid)
		dump += 'lcurve nser=%d ' % len(self.param['lc_energy_bands'])
		i = 1
		for eband in self.param['lc_energy_bands']:
			fenesel = 'ni%s_0mpu7_cl_%s.evt' % (self.obsid,eband) 
			dump += 'cfile%d="%s/%s" ' % (i,suboutdir,fenesel)
			i += 1 
		dump += 'window="-" ' 
		dump += 'dtnb=%d ' % self.param['lc_time_bin_sec']
		dump += 'nbint=%d ' % nbint
		dump += 'outfile="%s" ' % flcfile 
		dump += 'plotdnum=%s ' % len(self.param['lc_energy_bands'])
		dump += 'plot=yes plotdev="/xw" <<EOF\n'
		dump += 'lwid 5\n'
		dump += 'la ot %s\n' % title 
		dump += 'lab rotate\n'
		dump += 'lab pos y 2.8\n'
		i = 2
		for eband in self.param['lc_energy_bands']:
			emin,emax = self.param['lc_energy_bands'][eband]		
			dump += 'lab y%d %.1f-%.1f keV\n' % (i,emin,emax)
			dump += 'col %d on %d\n' % (i,i)
			i += 1 
		dump += 'hard ni%s_0mpu7_cl_ene.ps/cps\n' % self.obsid
		dump += 'quit\n'
		dump += 'EOF'
		print(dump)
		f.write(dump)
		f.close()

		# run the script.
		cmd  = 'chmod +x %s\n' % fcmd
		cmd += './%s' % fcmd
		print(cmd)
		os.system(cmd)

		cmd  = "ps2pdf ni%s_0mpu7_cl_ene.ps\n" % self.obsid
		cmd += "rm -f ni%s_0mpu7_cl_ene.ps\n" % self.obsid
		cmd += "mv ni%s_0mpu7_cl_ene.pdf %s/lc/\n" % (self.obsid,self.outdir)
		print(cmd)
		os.system(cmd)			

		self.dump_setup_to_yamlfile()

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

	def devide_to_gti(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		suboutdir = '%s/gti' % self.outdir
		cmd  = 'rm -rf %s;\n' % suboutdir
		cmd += 'mkdir -p %s\n' % suboutdir
		print(cmd);os.system(cmd)

		hdu = fits.open(self.clevt)
		number_of_gti = len(hdu['GTI'].data)

		for gtinum in range(number_of_gti):
			gti_dir = '%s/%003d' % (suboutdir,gtinum)
			cmd  = 'mkdir -p %s;\n' % gti_dir
			cmd += 'extract_gti_fits.py %s %d;\n' % (self.clevt,gtinum)
			gti_basename = 'ni%s_0mpu7_cl_gti%03d' % (self.obsid,gtinum)
			cmd += 'mv %s.{evt,gti} %s\n' % (gti_basename,gti_dir)
			print(cmd);os.system(cmd)

			nigti = NicerGTI(self,gtinum,gti_dir,gti_basename)
			if len(nigti.gti_hdu['EVENTS'].data) == 0:
				print("... No event. skip.")
				continue

			nigti.set_title()
			if nigti.nbint > 0: nigti.plot_lightcurve()
			nigti.run_nibackgen3C50()

class NicerGTI():
	def __init__(self,parent_nicerobsid,gtinum,gti_dir,gti_basename):
		self.parent_nicerobsid = parent_nicerobsid
		self.gtinum = gtinum; 
		self.gti_dir = gti_dir
		self.gti_basename = gti_basename
		print("\n[NicerGTI %d] A new NicerGTI is generated from a parent NicerObdID %s" % (
			self.gtinum,self.parent_nicerobsid.obsid))

		self.obsid = self.parent_nicerobsid.obsid
		self.param = self.parent_nicerobsid.param

		self.gti_fits = '%s/%s.gti' % (self.gti_dir,self.gti_basename)
		self.gti_clevt = '%s/%s.evt' % (self.gti_dir,self.gti_basename)
		self.gti_hdu = fits.open(self.gti_clevt)

	def set_title(self):
		print("\n[NicerGTI %s/%s] %s" % (self.gtinum,self.obsid,sys._getframe().f_code.co_name))

		self.tstart = float(self.gti_hdu['GTI'].data[0][0])
		self.tstop = float(self.gti_hdu['GTI'].data[0][1])
		self.exposure = self.tstop - self.tstart
		self.nbint = round(self.exposure/float(self.param['lc_time_bin_sec']))
		self.title = 'ObsID:%s GTI:%003d MJD:%.2f-%.2f (%.1f s)' % (
			self.obsid,self.gtinum,self.tstart,self.tstop,self.exposure)

	def plot_lightcurve(self):
		print("\n[NicerGTI %s/%s] %s" % (self.gtinum,self.obsid,sys._getframe().f_code.co_name))

		cmd = '\n'
		for eband in self.param['lc_energy_bands']:
			emin,emax = self.param['lc_energy_bands'][eband]
			fenesel = '%s/ni%s_0mpu7_cl_gti%03d_%s.evt' % (self.gti_dir,self.obsid,self.gtinum,eband)
			cmd += "fselect_energy.py %s %s %.1f %.1f\n" % (self.gti_clevt,fenesel,emin,emax)
		print(cmd);os.system(cmd)			

		fcmd = '%s/ni%s_0mpu7_cl_gti%03d_ene.sh' % (self.gti_dir,self.obsid,self.gtinum)
		flog = '%s/ni%s_0mpu7_cl_gti%03d_ene.log' % (self.gti_dir,self.obsid,self.gtinum)
		psfile = '%s/ni%s_0mpu7_cl_gti%03d_ene.ps' % (self.gti_dir,self.obsid,self.gtinum)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'

		flcfile = '%s/ni%s_0mpu7_cl_gti%03d_ene.flc' % (self.gti_dir,self.obsid,self.gtinum)
		dump += 'lcurve nser=%d ' % len(self.param['lc_energy_bands'])
		i = 1
		for eband in self.param['lc_energy_bands']:
			fenesel = '%s/ni%s_0mpu7_cl_gti%03d_%s.evt' % (self.gti_dir,self.obsid,self.gtinum,eband)
			dump += 'cfile%d="%s" ' % (i,fenesel)
			i += 1 
		dump += 'window="-" ' 
		dump += 'dtnb=%d ' % self.param['lc_time_bin_sec']
		dump += 'nbint=%d ' % self.nbint
		dump += 'outfile="%s" ' % flcfile 
		dump += 'plotdnum=%s ' % len(self.param['lc_energy_bands'])
		dump += 'plot=yes plotdev="/xw" <<EOF\n'
		dump += 'lwid 5\n'
		dump += 'la ot %s\n' % self.title 
		dump += 'lab rotate\n'
		dump += 'lab pos y 2.8\n'
		i = 2
		for eband in self.param['lc_energy_bands']:
			emin,emax = self.param['lc_energy_bands'][eband]		
			dump += 'lab y%d %.1f-%.1f keV\n' % (i,emin,emax)
			dump += 'col %d on %d\n' % (i,i)
			i += 1 
		dump += 'hard %s/cps\n' % psfile
		dump += 'quit\n'
		dump += 'EOF'
		print(dump)
		f.write(dump)
		f.close()	

		# run the script.
		cmd  = 'chmod +x %s\n' % fcmd
		cmd += './%s' % fcmd
		print(cmd)
		os.system(cmd)

		cmd  = "ps2pdf %s\n" % psfile
		cmd += "rm -f %s\n" % psfile
		cmd += "mv %s.pdf %s\n" % (os.path.splitext(os.path.basename(psfile))[0],self.gti_dir)
		print(cmd)
		os.system(cmd)

	def run_nibackgen3C50(self):
		print("\n[NicerGTI %s/%s] %s" % (self.gtinum,self.obsid,sys._getframe().f_code.co_name))

		self.gti_ufaevt = '%s/ni%s_0mpu7_ufa_gti%03d.evt' % (self.gti_dir,self.obsid,self.gtinum)

		cmd = 'rm -f xselect.log'
		print(cmd);os.system(cmd)

		cmd  = 'xselect<<EOF\n'
		cmd += 'xsel\n'
		cmd += 'read event %s ./\n' % self.parent_nicerobsid.ufaevt
		cmd += 'yes\n'
		cmd += 'filter time file\n'
		cmd += '%s\n' % self.gti_fits
		cmd += 'extract event\n'
		cmd += 'save event %s\n' % self.gti_ufaevt
		cmd += 'yes\n'
		cmd += 'exit\n'
		cmd += 'no\n'
		cmd += 'exit\n'	
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)	

		cmd = 'rm -f tmp_gti.txt xselect.log xsel_timefile.asc'
		print(cmd);os.system(cmd)

		# ============================================

		self.totspec_name = 'ni%s_gti%03d_3c50_tot' % (self.obsid,self.gtinum)
		self.bkgspec_name = 'ni%s_gti%03d_3c50_bkg' % (self.obsid,self.gtinum)
		self.totspec = '%s/%s.pi' % (self.gti_dir,self.totspec_name)
		self.bkgspec = '%s/%s.pi' % (self.gti_dir,self.bkgspec_name)
		# prepare a script for nicerl2 for each ObsID
		fcmd = '%s/nibackgen3C50_%s_gti%03d.sh' % (self.gti_dir,self.obsid,self.gtinum)
		flog = '%s/nibackgen3C50_%s_gti%03d.log' % (self.gti_dir,self.obsid,self.gtinum)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'
		dump += 'nibackgen3C50 '
		dump += 'rootdir=\'NONE\' '
		dump += 'obsid=\'NONE\' ' 
		dump += 'calevtdir=\'NONE\' ' 
		dump += 'bkgidxdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']
		dump += 'bkglibdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']	
		dump += 'gainepoch=\'%s\' ' % self.param['nibackgen3c50_gainepoch']
		dump += 'clfile=\'%s\' ' % self.gti_clevt		
		dump += 'ufafile=\'%s\' ' % self.gti_ufaevt
		dump += 'totspec=\'%s/%s\' ' % (self.gti_dir,self.totspec_name)
		dump += 'bkgspec=\'%s/%s\' ' % (self.gti_dir,self.bkgspec_name)
		dump += 'dtmin=%.1f dtmax=%.1f hbgcut=%.1f ' % (self.param['nibackgen3c50_dtmin'],self.param['nibackgen3c50_dtmax'],self.param['nibackgen3c50_hbgcut'])
		dump += '> %s 2>&1 ' % flog
		dump += '\n'
		f.write(dump)
		f.close()	

		try:
			# run the script.
			cmd  = 'chmod +x %s\n' % fcmd
			cmd += './%s' % fcmd
			print(cmd);os.system(cmd)	
		except:
			f = open(flog,'w')
			message = "nibackgen3c50 failed at gti %d (e.g., 'Error of no events in new cleaned event file')" % self.gtinum
			print(message)
			f.write(message)
			f.close()
			return -1

		if not (os.path.exists(self.totspec) and os.path.exists(self.bkgspec)):
			f = open(flog,'w')
			message = "nibackgen3c50 failed at gti %d (e.g., 'Error of no events in new cleaned event file')" % self.gtinum
			print(message)
			f.write(message)
			f.close()
			return -1			

		cmd  = 'fparkey %s %s RESPFILE\n' % (self.param['xspec_rmf'],self.totspec)
		cmd += 'fparkey %s %s ANCRFILE\n' % (self.param['xspec_arf'],self.totspec)
		#cmd += 'fparkey %s %s BACKFILE\n' % (bkgspec,totspec)
		cmd += 'fparkey %s %s RESPFILE\n' % (self.param['xspec_rmf'],self.bkgspec)
		cmd += 'fparkey %s %s ANCRFILE\n' % (self.param['xspec_arf'],self.bkgspec)
		print(cmd);os.system(cmd)

		hdu = fits.open(self.totspec)
		keyword_dateobs = hdu[0].header['DATE-OBS']
		keyword_object = hdu[0].header['OBJECT']		
		keyword_exposure = hdu[0].header['EXPOSURE']			
		title = '%s ObsID:%s GTI:%03d (%s, %.1f sec)' % (keyword_object, self.obsid, self.gtinum, keyword_dateobs,keyword_exposure)

		fspec = "ni%s_tot_spec" % self.obsid
		cmd  = "xspec <<EOF\n"
		cmd += "data 1 %s\n" % (self.totspec)
		cmd += "back 1 %s\n" % (self.bkgspec)
		cmd += "data 2 %s\n" % (self.bkgspec)
		cmd += "data 3 %s\n" % (self.totspec)	
		cmd += "back 3 none\n"
		cmd += "setplot energy\n"
		cmd += "setplot rebin %d %d\n" % (self.param['xspec_rebin_sigma'],self.param['xspec_rebin_maxnum'])
		cmd += "ignore 1-3:**-%.1f,%.1f-**\n" % (self.param['xspec_emin'],self.param['xspec_emax'])
		cmd += "iplot %s\n" % self.param['xspec_plot_type']
		cmd += "time off\n"
		cmd += "lwid 5\n"
		cmd += "lwid 5 on 1..3\n"
		cmd += "la t %s\n" % title
		cmd += "r x %.1f %.1f\n" % (self.param['xspec_emin'],self.param['xspec_emax'])
		#cmd += "r y 1e-4 100.0\n"
		cmd += "r y %.1f %.1f\n" % (self.param['xspec_ymin'],self.param['xspec_ymax'])
		cmd += "hard %s.ps/cps\n" % fspec
		cmd += "exit\n"
		cmd += "exit\n"
		cmd += "EOF\n"
		print(cmd)
		os.system(cmd)

		cmd  = "ps2pdf %s.ps\n" % fspec
		cmd += "rm -f %s.ps\n" % fspec
		cmd += "mv %s.pdf %s/\n" % (fspec,self.gti_dir)
		print(cmd)
		os.system(cmd)


class NicerElf():
	def __init__(self,setup_yamlfile,obsid_lstfile):
		print("\n[NicerElf] A new NicerElf object is generated.")

		print('...setup_yamlfile: {}'.format(setup_yamlfile))
		print('...obsid_lstfile: {}'.format(obsid_lstfile))
		self.param = yaml.load(open(setup_yamlfile),Loader=yaml.SafeLoader)
		self.obsid_lstfile = obsid_lstfile

		print('...param: %s' % self.param)

		print('...seting nicerobs_lst')
		self.nicerobs_lst = []
		for line in open(self.obsid_lstfile):
			obsid = line.split()[0]
			if obsid in ['#','%']:
				print('...comment out and skip the line: %s' % line)
				continue 
			elif obsid == 'exit':
				break
			elif len(glob.glob('%s/%s' % (self.param["input_data_directory"],obsid))) == 0:
				print('...no directory and skip obsid: %s' % obsid)
				continue
			else:
				indir = glob.glob('%s/%s' % (self.param["input_data_directory"],obsid))[0]
				outdir = '%s/%s' % (self.param['output_directory'],obsid)
				print('...set obsid:%s' % indir)
				self.nicerobs_lst.append(NicerObsID(
					obsid=obsid,indir=indir,outdir=outdir,param=self.param))

	def make_directory(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_make_directory']:
			print("skip ... since flag_make_directory is %s" % self.param['flag_make_directory']) 
			return 0		
		for niobsid in self.nicerobs_lst:
			niobsid.make_directory()

	def run_nicerl2(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nicerl2']:
			print("skip ... since flag_run_nicerl2 is %s" % self.param['flag_run_nicerl2']) 
			return 0			
		for niobsid in self.nicerobs_lst:
			niobsid.run_nicerl2()

	def run_nibackgen3C50(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nibackgen3C50']:
			print("skip ... since flag_run_nibackgen3C50 is %s" % self.param['flag_run_nibackgen3C50']) 
			return 0			
		for niobsid in self.nicerobs_lst:
			niobsid.run_nibackgen3C50()

	def plot_lightcurve(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_plot_lightcurve']:
			print("skip ... since flag_plot_lightcurve is %s" % self.param['flag_plot_lightcurve']) 
			return 0			
		for niobsid in self.nicerobs_lst:
			niobsid.plot_lightcurve()

	def run_barycorr(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_barycorr']:
			print("skip ... since flag_run_barycorr is %s" % self.param['flag_run_barycorr']) 
			return 0			
		for niobsid in self.nicerobs_lst:
			niobsid.run_barycorr()

	def devide_to_gti(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_devide_to_gti']:
			print("skip ... since flag_devide_to_gti is %s" % self.param['flag_devide_to_gti']) 
			return 0		
		for niobsid in self.nicerobs_lst:
			niobsid.devide_to_gti()

"""
def record_results(obsid):
	print("record_results obsid=%s" % obsid)

	if not os.path.exists('%s/%s' % (outdir,obsid)):
		print("ERROR no directory for %s." % obsid)
		return -1

	cmd  = 'rm -rf %s/%s/param;\n' % (outdir,obsid)
	cmd += 'mkdir -p %s/%s/param\n' % (outdir,obsid)
	print(cmd);os.system(cmd)

	totspec_name = 'ni%s_3c50_tot' % obsid
	bkgspec_name = 'ni%s_3c50_bkg' % obsid	
	totspec = '%s/%s/spec/%s.pi' % (outdir,obsid,totspec_name)
	bkgspec = '%s/%s/spec/%s.pi' % (outdir,obsid,bkgspec_name)

	param_dict = {}
	param_dict["ratechk_total"] = {}
	for eband in param['ratechk_energy_bands']:
		emin,emax = param['ratechk_energy_bands'][eband]
		cmd  = 'rm -f tmp.log;\n'
		cmd += "extract_xspec_rate.py "
		cmd += "%s " % totspec
		cmd += "%s " % param['xspec_rmf']
		cmd += "%s " % param['xspec_arf']	
		cmd += "%s " % emin
		cmd += "%s " % emax
		cmd += "--keyword R%s " % eband	 	
		cmd += "> tmp.log\n"
		print(cmd)
		os.system(cmd)

		#with open('tmp.log', mode='rt', encoding='utf-8') as f:
		#	read_data = list(f)
		#src_rate  = float(read_data[-2].split(":")[-1])
		##src_error = float(read_data[-1].split(":")[-1])	
		#print("src",src_rate,src_error)

		cmd  = 'rm -f tmp.log;\n'
		cmd += "extract_xspec_rate.py "
		cmd += "%s " % bkgspec
		cmd += "%s " % param['xspec_rmf']
		cmd += "%s " % param['xspec_arf']	
		cmd += "%s " % emin		
		cmd += "%s " % emax	
		cmd += "--keyword R%s " % eband	 
		cmd += "> tmp.log\n"
		print(cmd)
		os.system(cmd)

		#with open('tmp.log', mode='rt', encoding='utf-8') as f:
		#	read_data = list(f)
		#bkg_rate  = float(read_data[-2].split(":")[-1])
		#bkg_error = float(read_data[-1].split(":")[-1])	
		#print("bkg",bkg_rate,bkg_error)	

		#sub_rate = float('%.6f' % (src_rate - bkg_rate))
		#sub_error = float('%.6f' % (np.sqrt(src_error*src_error+bkg_error*bkg_error)))
		#print("sub",sub_rate,sub_error)

		#param_dict["ratechk_total"][eband] = [src_rate,src_error,bkg_rate,bkg_error,sub_rate,sub_error]

		cmd = 'rm -f tmp.log'
		print(cmd);os.system(cmd)

	fparam = '%s/%s/param/ni%s.yaml' % (outdir,obsid,obsid)
	with open(fparam, 'w') as file:
		documents = yaml.dump(param_dict, file)	
"""		