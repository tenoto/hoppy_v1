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


		self.error_log = '%s/error_log.txt' % self.param['output_directory']

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
		self.spec_pdf = '%s/%s.pdf' % (suboutdir,fspec)
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

		self.lc_pdf = '%s/lc/ni%s_0mpu7_cl_ene.pdf' % (self.outdir,self.obsid)

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

	def devide_to_segment(self):
		print("\n[NicerObsID %s] %s" % (self.obsid,sys._getframe().f_code.co_name))

		suboutdir = '%s/segment' % self.outdir
		cmd  = 'rm -rf %s;\n' % suboutdir
		cmd += 'mkdir -p %s\n' % suboutdir
		print(cmd);os.system(cmd)

		print(self.clevt)
		cmd = 'rm -f tmp_gti.txt; group_gti.py %s tmp_gti.txt --gti_separation_threshold %.2f' % (self.clevt,self.param['gti_separation_threshold'])
		print(cmd);os.system(cmd)

		segment_gti = '%s/ni%s_segment.gti' % (suboutdir,self.obsid)
		cmd = 'fconv_txt2gti.py -i tmp_gti.txt -o %s; rm -f tmp_gti.txt' % segment_gti 
		print(cmd);os.system(cmd)

		hdu = fits.open(segment_gti)
		number_of_segment = len(hdu['GTI'].data)

		self.nigsegment_list = []
		for i in range(number_of_segment):
			segment_num = i + 1 
			segment_dir = '%s/ep%003d' % (suboutdir,segment_num)
			print("----- segment_dir = %s" % segment_dir)
			cmd  = 'mkdir -p %s;\n' % segment_dir
			cmd += 'cp %s %s/tmp_sel.gti;\n' % (segment_gti,segment_dir)
			if segment_num == 1:
				if number_of_segment == 2:
					#cmd += 'ftdelrow %s/tmp_sel.gti+1 none 2-2 confirm=YES;\n' % (segment_dir)
					cmd += 'fdelrow %s/tmp_sel.gti+1 2 1 Y Y;\n' % (segment_dir)
				else:
					cmd += 'ftdelrow %s/tmp_sel.gti+1 none %d- confirm=YES;\n' % (segment_dir,segment_num+1)					
			elif segment_num == number_of_segment:
				if number_of_segment == 2:
					cmd += 'ftdelrow %s/tmp_sel.gti+1 none 1 confirm=YES;\n' % (segment_dir)
				else:
					cmd += 'ftdelrow %s/tmp_sel.gti+1 none 1-%d confirm=YES;\n' % (segment_dir,segment_num-1)	
			else:
				cmd += 'ftdelrow %s/tmp_sel.gti+1 none 1-%d,%d- confirm=YES;\n' % (segment_dir,segment_num-1,segment_num+1)
			print(cmd);os.system(cmd)

			segment_basename = 'ni%s_0mpu7_cl_ep%03d' % (self.obsid,segment_num)
			segment_fitsfile = '%s.evt' % segment_basename
			cmd  = 'xselect_gtifilter.py '
			cmd += '-i %s ' % self.clevt
			cmd += '-g %s/tmp_sel.gti ' % segment_dir
			cmd += '-o %s/%s;\n' % (segment_dir,segment_fitsfile)
			print(cmd);os.system(cmd)

			nisegment = NicerSegment(self,segment_num,segment_dir,segment_basename)
			if len(nisegment.segment_clhdu['EVENTS'].data) == 0:
				print("... No event. skip.")
				continue

			nisegment.set_title()
			if nisegment.nbint > 0: 
				nisegment.plot_lightcurve()
			nisegment.run_nibackgen3C50()
			self.nigsegment_list.append(nisegment)

		self.df_summary = pd.DataFrame(columns=['ObsID','segment','lc','spec'])
		#for i in range(number_of_segment):
		for nisegment in self.nigsegment_list:
			lc_link = '<a href="./%s">pdf</a>' % nisegment.lc_pdf.replace(suboutdir,'')
			spec_link = '<a href="./%s">pdf</a>' % nisegment.spec_pdf.replace(suboutdir,'')			
			self.df_summary = self.df_summary.append({'ObsID':self.obsid,
				'segment':nisegment.segment_num,
				'lc':lc_link,'spec':spec_link}, ignore_index=True)

		self.show_dataframe_summay()
		self.write_dataframe_summary()

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
				f = open(self.error_log,'w')
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


class NicerSegment():
	def __init__(self,parent_nicerobsid,segment_num,segment_dir,segment_basename):
		self.parent_nicerobsid = parent_nicerobsid
		self.segment_num = segment_num; 
		self.segment_dir = segment_dir
		self.segment_basename = segment_basename
		print("\n[NicerSegment %d] A new NicerSegment is generated from a parent NicerObdID %s" % (
			self.segment_num,self.parent_nicerobsid.obsid))

		self.obsid = self.parent_nicerobsid.obsid
		self.param = self.parent_nicerobsid.param

		self.segment_fits = '%s/%s.gti' % (self.segment_dir,self.segment_basename)
		self.segment_clevt = '%s/%s.evt' % (self.segment_dir,self.segment_basename)
		self.segment_clhdu = fits.open(self.segment_clevt)

		self.segment_ufaevt = '%s/ni%s_0mpu7_ufa_ep%03d.evt' % (self.segment_dir,self.obsid,self.segment_num)
		self.totspec_name = 'ni%s_ep%03d_3c50_tot' % (self.obsid,self.segment_num)
		self.bkgspec_name = 'ni%s_ep%03d_3c50_bkg' % (self.obsid,self.segment_num)
		self.totspec = '%s/%s.pi' % (self.segment_dir,self.totspec_name)
		self.bkgspec = '%s/%s.pi' % (self.segment_dir,self.bkgspec_name)

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

	def plot_lightcurve(self):
		print("\n[NicerSegment %s/%s] %s" % (self.segment_num,self.obsid,sys._getframe().f_code.co_name))

		cmd = '\n'
		for eband in self.param['lc_energy_bands']:
			emin,emax = self.param['lc_energy_bands'][eband]
			fenesel = '%s/ni%s_0mpu7_cl_ep%03d_%s.evt' % (self.segment_dir,self.obsid,self.segment_num,eband)
			cmd += "fselect_energy.py %s %s %.1f %.1f\n" % (self.segment_clevt,fenesel,emin,emax)
		print(cmd);os.system(cmd)			

		fcmd = '%s/ni%s_0mpu7_cl_ep%03d_ene.sh' % (self.segment_dir,self.obsid,self.segment_num)
		flog = '%s/ni%s_0mpu7_cl_ep%03d_ene.log' % (self.segment_dir,self.obsid,self.segment_num)
		psfile = '%s/ni%s_0mpu7_cl_ep%03d_ene.ps' % (self.segment_dir,self.obsid,self.segment_num)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'

		flcfile = '%s/ni%s_0mpu7_cl_ep%03d_ene.flc' % (self.segment_dir,self.obsid,self.segment_num)
		dump += 'lcurve nser=%d ' % len(self.param['lc_energy_bands'])
		i = 1
		for eband in self.param['lc_energy_bands']:
			fenesel = '%s/ni%s_0mpu7_cl_ep%03d_%s.evt' % (self.segment_dir,self.obsid,self.segment_num,eband)
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
		cmd += "mv %s.pdf %s\n" % (os.path.splitext(os.path.basename(psfile))[0],self.segment_dir)
		print(cmd)
		os.system(cmd)

		self.lc_pdf = '%s/%s.pdf' % (self.segment_dir,os.path.splitext(os.path.basename(psfile))[0])

	def run_nibackgen3C50(self):
		print("\n[NicerSegment %s/%s] %s" % (self.segment_num,self.obsid,sys._getframe().f_code.co_name))

		cmd = 'rm -f xselect.log'
		print(cmd);os.system(cmd)

		cmd  = 'xselect<<EOF\n'
		cmd += 'xsel\n'
		cmd += 'read event %s ./\n' % self.parent_nicerobsid.ufaevt
		cmd += 'yes\n'
		cmd += 'filter time file\n'
		cmd += '%s\n' % self.segment_fits
		cmd += 'extract event\n'
		cmd += 'save event %s\n' % self.segment_ufaevt
		cmd += 'yes\n'
		cmd += 'exit\n'
		cmd += 'no\n'
		cmd += 'exit\n'	
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)	

		cmd = 'rm -f tmp_gti.txt xselect.log xsel_timefile.asc'
		print(cmd);os.system(cmd)

		# ============================================


		# prepare a script for nicerl2 for each ObsID
		fcmd = '%s/nibackgen3C50_%s_ep%03d.sh' % (self.segment_dir,self.obsid,self.segment_num)
		flog = '%s/nibackgen3C50_%s_ep%03d.log' % (self.segment_dir,self.obsid,self.segment_num)
		f = open(fcmd,'w')
		dump  = '#!/bin/sh -f\n'
		dump += 'nibackgen3C50 '
		dump += 'rootdir=\'NONE\' '
		dump += 'obsid=\'NONE\' ' 
		dump += 'calevtdir=\'NONE\' ' 
		dump += 'bkgidxdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']
		dump += 'bkglibdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']	
		dump += 'gainepoch=\'%s\' ' % self.param['nibackgen3c50_gainepoch']
		dump += 'clfile=\'%s\' ' % self.segment_clevt		
		dump += 'ufafile=\'%s\' ' % self.segment_ufaevt
		dump += 'totspec=\'%s/%s\' ' % (self.segment_dir,self.totspec_name)
		dump += 'bkgspec=\'%s/%s\' ' % (self.segment_dir,self.bkgspec_name)
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
			message = "nibackgen3c50 failed at ep %d (e.g., 'Error of no events in new cleaned event file')" % self.segment_num
			print(message)
			f.write(message)
			f.close()
			return -1

		if not (os.path.exists(self.totspec) and os.path.exists(self.bkgspec)):
			f = open(flog,'w')
			message = "nibackgen3c50 failed at ep %d (e.g., 'Error of no events in new cleaned event file')" % self.segment_num
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
		title = '%s ObsID:%s GTI:%03d (%s, %.1f sec)' % (keyword_object, self.obsid, self.segment_num, keyword_dateobs,keyword_exposure)

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
		cmd += "mv %s.pdf %s/\n" % (fspec,self.segment_dir)
		print(cmd)
		os.system(cmd)

		self.spec_pdf = '%s/%s.pdf' % (self.segment_dir,fspec)

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

		#self.fit_pdf = '%s/fit/%s.pdf' % (self.segment_dir,fspec)


class NicerElf():
	def __init__(self,setup_yamlfile,obsid_lstfile):
		print("\n[NicerElf] A new NicerElf object is generated.")

		print('...setup_yamlfile: {}'.format(setup_yamlfile))
		print('...obsid_lstfile: {}'.format(obsid_lstfile))
		self.param = yaml.load(open(setup_yamlfile),Loader=yaml.SafeLoader)
		self.obsid_lstfile = obsid_lstfile

		print('...param: %s' % self.param)

		self.error_log = '%s/error_log.txt' % self.param['output_directory']
		self.df_summary = pd.DataFrame(columns=['ObsID'])

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
					break 

			if obsid in ['#','%']:
				print('...comment out and skip the line: %s' % line)
				continue 
			elif obsid == 'exit' or obsid == 'quit' or obsid == 'exit()' or obsid == 'quit()' :
				break
			elif len(glob.glob('%s/%s' % (self.param["input_data_directory"],obsid))) == 0:
				print('...no directory and skip obsid: %s' % obsid)
				continue					
			elif not flag_have_clevt:
				print('...no clenaed event (gpg?) and skip obsid: %s' % obsid)
				continue					
			else:
				indir = glob.glob('%s/%s' % (self.param["input_data_directory"],obsid))[0]
				outdir = '%s/%s' % (self.param['output_directory'],obsid)
				print('...set obsid:%s' % indir)
				self.nicerobs_lst.append(NicerObsID(
					obsid=obsid,indir=indir,outdir=outdir,param=self.param))

		for niobsid in self.nicerobs_lst:
			self.df_summary = self.df_summary.append({'ObsID':niobsid.obsid}, ignore_index=True)

	def show_dataframe_summay(self):
		print("******* summary ******")
		print(self.df_summary)
		print("**********************")		

	def write_dataframe_summary(self):
		self.summary_csv = '%s/nipipeline_summary.csv' % self.param['output_directory']
		self.summary_html = '%s/nipipeline_summary.html' % self.param['output_directory']		
		self.df_summary.to_csv(self.summary_csv)

		raw_html = self.df_summary.to_html(open(self.summary_html,'w'))
		cmd = 'sed -e "s/&lt;/</g" %s > tmp.html' % self.summary_html
		print(cmd);os.system(cmd)
		cmd = 'sed -e "s/&gt;/>/g" tmp.html > tmp2.html' 
		print(cmd);os.system(cmd)		
		cmd = 'mv tmp2.html %s; rm -f tmp.html tmp2.html' % self.summary_html
		print(cmd);os.system(cmd)		

	def make_directory(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_make_directory']:
			print("skip ... since flag_make_directory is %s" % self.param['flag_make_directory']) 
			return 0		
		dict_append = {}
		for niobsid in self.nicerobs_lst:
			niobsid.make_directory()
			link = '<a href="./%s">dir</a>' % niobsid.outdir.replace(self.param['output_directory'],'')
			dict_append[link] = niobsid.obsid
		self.df_summary['Directory'] = dict_append
		self.write_dataframe_summary()

	def run_nicerl2(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nicerl2']:
			print("skip ... since flag_run_nicerl2 is %s" % self.param['flag_run_nicerl2']) 
			return 0
		dict_append = {}						
		for niobsid in self.nicerobs_lst:
			niobsid.run_nicerl2()
			link = '<a href="./%s">log</a>' % niobsid.flog_nicerl2.replace(self.param['output_directory'],'')
			dict_append[link] = niobsid.obsid
		self.df_summary['nicerl2'] = dict_append
		self.write_dataframe_summary()

	def run_nibackgen3C50(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_nibackgen3C50']:
			print("skip ... since flag_run_nibackgen3C50 is %s" % self.param['flag_run_nibackgen3C50']) 
			return 0	
		dict_append = {}											
		for niobsid in self.nicerobs_lst:
			niobsid.run_nibackgen3C50()
			link = '<a href="./%s">pdf</a>' % niobsid.spec_pdf.replace(self.param['output_directory'],'')
			dict_append[link] = niobsid.obsid
		self.df_summary['spec'] = dict_append
		self.write_dataframe_summary()

	def plot_lightcurve(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_plot_lightcurve']:
			print("skip ... since flag_plot_lightcurve is %s" % self.param['flag_plot_lightcurve']) 
			return 0	
		dict_append = {}														
		for niobsid in self.nicerobs_lst:
			niobsid.plot_lightcurve()
			link = '<a href="./%s">pdf</a>' % niobsid.lc_pdf.replace(self.param['output_directory'],'')
			dict_append[link] = niobsid.obsid
		self.df_summary['lc'] = dict_append
		self.write_dataframe_summary()

	def run_barycorr(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_run_barycorr']:
			print("skip ... since flag_run_barycorr is %s" % self.param['flag_run_barycorr']) 
			return 0			
		for niobsid in self.nicerobs_lst:
			niobsid.run_barycorr()

	def devide_to_segment(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_devide_to_segment']:
			print("skip ... since flag_devide_to_segment is %s" % self.param['flag_devide_to_segment']) 
			return 0		
		dict_append = {}				
		for niobsid in self.nicerobs_lst:
			niobsid.devide_to_segment()
			link = '<a href="./%s">html</a>' % niobsid.summary_html.replace(self.param['output_directory'],'')
			dict_append[link] = niobsid.obsid
		self.df_summary['segment'] = dict_append
		self.write_dataframe_summary()

	def fit_of_segment(self):
		print("\n[NicerElf] %s" % sys._getframe().f_code.co_name)
		if not self.param['flag_fit_of_segment']:
			print("skip ... since flag_fit_of_segment is %s" % self.param['flag_fit_of_segment']) 
			return 0		
		for niobsid in self.nicerobs_lst:
			niobsid.fit_of_segment()
