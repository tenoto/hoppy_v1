# -*- coding: utf-8 -*-

import os 
import glob
import yaml 
import pandas as pd

from astropy.io import fits
from astropy.time import Time

class NicerElf():
	def __init__(self,setup_yamlfile,obsid_lstfile):
		print("... NicerElf object is generated")
		self.param = yaml.load(open(setup_yamlfile),Loader=yaml.SafeLoader)
		self.obsid_lstfile = obsid_lstfile

		self.obsid_lst = []
		for line in open(self.obsid_lstfile):
			obsid_text = line.split()[0]
			if obsid_text in ['#','%']:
				print('skip')
				continue 
			for indir in glob.glob('%s/%s' % (self.param["input_data_directory"],obsid_text)):
				obsid = os.path.basename(indir)
				print(obsid,indir)
				self.obsid_lst.append(obsid)

		print(self.param)
		print(self.obsid_lst)

	def make_directory(self):
		print("...make_directory")
		for obsid in self.obsid_lst:
			print("make_directory obsid=%s" % obsid)
			indir = glob.glob('%s/%s' % (self.param['input_data_directory'],obsid))[0]
			outdir = '%s/%s' % (self.param['output_directory'],obsid)

			cmd  = 'rm -rf %s;\n' % (outdir)
			cmd += 'mkdir -p %s/xti;\n' % (outdir)
			cmd += 'ln -s %s/{auxil,log} %s/;\n' % (indir,outdir)
			cmd += 'ln -s %s/xti/{event_uf,hk} %s/xti/;\n' % (indir,outdir)	
			print(cmd);os.system(cmd)

	def run_nicerl2(self):
		print("...run_nicerl2")		
		fversion = os.popen('fversion').read().rstrip()		
		for obsid in self.obsid_lst:
			print("run_nicerl2 obsid=%s" % obsid)
			outdir = '%s/%s' % (self.param['output_directory'],obsid)
			
			# prepare a script for nicerl2 for each ObsID
			fcmd = '%s/nicerl2_%s.sh' % (outdir,obsid)
			flog = '%s/nicerl2_%s.log' % (outdir,obsid)
			f = open(fcmd,'w')
			dump  = '#!/bin/sh -f\n'
			dump += 'nicerl2 indir=%s ' % (outdir)
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

	def run_nibackgen3C50(self):
		print("...run_nibackgen3C50")		
		fversion = os.popen('fversion').read().rstrip()		

		for obsid in self.obsid_lst:
			print("run_nibackgen3C50 obsid=%s" % obsid)
			outdir = '%s/%s' % (self.param['output_directory'],obsid)

			cmd  = 'rm -rf %s/spec;\n' % (outdir)
			cmd += 'mkdir -p %s/spec\n' % (outdir)
			print(cmd);os.system(cmd)

			clevt = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (outdir,obsid)
			hdu = fits.open(clevt)
			if len(hdu['EVENTS'].data) == 0:
				flog = '%s/spec/nibackgen3C50_%s.log' % (outdir,obsid)
				f = open(flog,'w')
				message = 'no cleaned event in the original file %s' % clevt
				print(message)
				f.write(message)
				f.close()
				continue 

			totspec_name = 'ni%s_3c50_tot' % obsid
			bkgspec_name = 'ni%s_3c50_bkg' % obsid	

			# prepare a script for nicerl2 for each ObsID
			fcmd = '%s/spec/nibackgen3C50_%s.sh' % (outdir,obsid)
			flog = '%s/spec/nibackgen3C50_%s.log' % (outdir,obsid)
			f = open(fcmd,'w')
			dump  = '#!/bin/sh -f\n'
			dump += 'nibackgen3C50 ' 
			dump += 'rootdir=\'%s\' ' % self.param['output_directory']
			dump += 'obsid=\'%s\' ' % obsid
			dump += 'bkgidxdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']
			dump += 'bkglibdir=\'%s\' ' % self.param['nibackgen3c50_bkgidxdir']	
			dump += 'gainepoch=\'%s\' ' % self.param['nibackgen3c50_gainepoch']
			dump += 'totspec=\'./%s\' ' % (totspec_name)
			dump += 'bkgspec=\'./%s\' ' % (bkgspec_name)
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

			cmd = 'mv %s*.pi %s*.pi %s/spec' % (totspec_name,bkgspec_name,outdir)
			print(cmd);os.system(cmd)


			totspec = '%s/spec/%s.pi' % (outdir,totspec_name)
			bkgspec = '%s/spec/%s.pi' % (outdir,bkgspec_name)

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
			title = '%s ObsID:%s (%s, %.1f sec)' % (keyword_object, obsid,keyword_dateobs,keyword_exposure)

			fspec = "ni%s_tot_spec" % obsid
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
			cmd += "mv %s.pdf %s/spec/\n" % (fspec,outdir)
			print(cmd)
			os.system(cmd)

"""
def plot_lightcurve(obsid):
	print("plot_lightcurve obsid=%s" % obsid)

	if not os.path.exists('%s/%s' % (outdir,obsid)):
		print("ERROR no directory for %s." % obsid)
		return -1

	cmd  = 'rm -rf %s/%s/lc;\n' % (outdir,obsid)
	cmd += 'mkdir -p %s/%s/lc\n' % (outdir,obsid)
	print(cmd);os.system(cmd)

	for eband in param['lc_energy_bands']:
		emin,emax = param['lc_energy_bands'][eband]
		fenesel = 'ni%s_0mpu7_cl_%s.evt' % (obsid,eband)
		cmd += "fselect_energy.py %s/%s/xti/event_cl/ni%s_0mpu7_cl.evt %s %.1f %.1f;\n" % (outdir,obsid,obsid,fenesel,emin,emax)
		cmd += "mv %s %s/%s/lc;\n" % (fenesel,outdir,obsid)
	print(cmd)
	os.system(cmd)

	hdu = fits.open('%s/%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (outdir,obsid,obsid))
	t = Time([hdu[1].header['DATE-OBS'],hdu[1].header['DATE-END']],format='isot', scale='utc')
	key_date_obs = hdu[1].header['DATE-OBS']      
	key_date_end = hdu[1].header['DATE-END']      
	key_exposure = hdu[1].header['EXPOSURE']
	key_target = hdu[1].header['OBJECT']
	key_mjd_start = t.mjd[0]
	key_mjd_stop = t.mjd[1]       
	key_onsource = (key_mjd_stop-key_mjd_start)*24.0*60.0*60.0
	nbint = round(float(key_onsource) / float(param['lc_time_bin_sec']))
	title = 'ObsID:%s %s MJD:%.2f (exposure %.1f s)' % (obsid,key_date_obs,key_mjd_start,key_exposure)

	fcmd = '%s/%s/lc/ni%s_0mpu7_cl_ene.sh' % (outdir,obsid,obsid)
	flog = '%s/%s/lc/ni%s_0mpu7_cl_ene.log' % (outdir,obsid,obsid)
	f = open(fcmd,'w')
	dump  = '#!/bin/sh -f\n'

	flcfile = '%s/%s/lc/ni%s_0mpu7_cl_ene.flc' % (outdir,obsid,obsid)
	dump += 'lcurve nser=%d ' % len(param['lc_energy_bands'])
	i = 1
	for eband in param['lc_energy_bands']:
		fenesel = 'ni%s_0mpu7_cl_%s.evt' % (obsid,eband) 
		dump += 'cfile%d="%s/%s/lc/%s" ' % (i,outdir,obsid,fenesel)
		i += 1 
	dump += 'window="-" ' 
	dump += 'dtnb=%d ' % param['lc_time_bin_sec']
	dump += 'nbint=%d ' % nbint
	dump += 'outfile="%s" ' % flcfile 
	dump += 'plotdnum=%s ' % len(param['lc_energy_bands'])
	dump += 'plot=yes plotdev="/xw" <<EOF\n'
	dump += 'lwid 5\n'
	dump += 'la ot %s\n' % title 
	dump += 'lab rotate\n'
	dump += 'lab pos y 2.8\n'
	i = 2
	for eband in param['lc_energy_bands']:
		emin,emax = param['lc_energy_bands'][eband]		
		dump += 'lab y%d %.1f-%.1f keV\n' % (i,emin,emax)
		dump += 'col %d on %d\n' % (i,i)
		i += 1 
	dump += 'hard ni%s_0mpu7_cl_ene.ps/cps\n' % obsid
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

	cmd  = "ps2pdf ni%s_0mpu7_cl_ene.ps\n" % obsid
	cmd += "rm -f ni%s_0mpu7_cl_ene.ps\n" % obsid
	cmd += "mv ni%s_0mpu7_cl_ene.pdf %s/%s/lc/\n" % (obsid,outdir,obsid)
	print(cmd)
	os.system(cmd)


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