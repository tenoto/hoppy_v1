# -*- coding: utf-8 -*-

import os 
import sys 
import yaml
import glob 
import pandas as pd
import numpy as np
import astropy.io.fits as fits 
import astropy.time as asttime 

def string_to_list(string):
	mainlist = []
	for col in string.split(','):
		sublist = []
		for i in col.split('-'):
			sublist.append(float(i))
		mainlist.append(sublist)
	return mainlist

def grep(logfile,keyword,colnum,dtype=float):
	out_word_list = []
	for line in open(logfile):
		if keyword in line:
			cols = line.split()
			if dtype == float:
				out_word_list.append(float(cols[colnum]))
			elif dtype == int:
				out_word_list.append(int(cols[colnum]))				
	return out_word_list

class XspecPha():
	def __init__(self,phafile,outdir='out',
		backgrnd=None,rmffile=None,arffile=None,modelxcm=None,
		binminsig=5,binmaxbin=50,fitemin=0.4,fitemax=10.0,
		plotxmin=0.4,plotxmax=10.0,plotymin=1e-5,plotymax=1e+4,
		plotymineeuf=1e-5,plotymaxeeuf=10,ploty2min=-10.0,ploty2max=10.0,
		ratebands=[[0.4,6.0],[1.0,10.0]],
		fluxbands=[[0.4,6.0],[1.0,10.0]],
		parerrnum=[1,2,5]):
		print("\n[XspecPha] %s" % (sys._getframe().f_code.co_name))

		self.phafile = phafile
		self.outdir = outdir 
		self.backgrnd = backgrnd
		self.rmffile = rmffile
		self.arffile = arffile
		self.modelxcm = modelxcm
		self.binminsig = binminsig
		self.binmaxbin = binmaxbin
		self.fitemin = fitemin
		self.fitemax = fitemax 
		self.plotxmin = plotxmin
		self.plotxmax = plotxmax
		self.plotymin = plotymin
		self.plotymax = plotymax 
		self.plotymineeuf = plotymineeuf
		self.plotymaxeeuf = plotymaxeeuf		
		self.ploty2min = ploty2min
		self.ploty2max = ploty2max

		self.flag_fit_complete = False

		if type(ratebands) == list:
			self.ratebands = ratebands
		elif type(ratebands) == str:
			self.ratebands = string_to_list(ratebands)

		if type(fluxbands) == list:
			self.fluxbands = fluxbands
		elif type(fluxbands) == str:
			self.fluxbands = string_to_list(fluxbands)		

		if type(parerrnum) == list:
			self.parerrnum = parerrnum
		elif type(fluxbands) == str:
			self.parerrnum = [int(i) for i in parerrnum.split(",")]

		if not os.path.exists(self.phafile):
			sys.stderr.write('phafile %s does not exist.' % self.phafile)
			exit()
		if os.path.exists(self.outdir):
			sys.stderr.write('output directory %s has already existed.' % self.outdir)
			quit()
		if not os.path.exists(self.outdir):
			cmd = 'mkdir -p %s' % self.outdir
			print(cmd);os.system(cmd)

		cmd = 'cp %s %s' % (self.phafile,self.outdir)
		print(cmd);os.system(cmd)
		self.phafile = '%s/%s' % (self.outdir,os.path.basename(self.phafile))

		if self.rmffile != None:
			cmd = "fparkey %s %s+1 RESPFILE" % (self.rmffile,self.phafile)
			print(cmd);os.system(cmd)
		if self.arffile != None:
			cmd = "fparkey %s %s+1 ANCRFILE" % (self.arffile,self.phafile)
			print(cmd);os.system(cmd)		

		if self.backgrnd != None:
			cmd = 'cp %s %s' % (self.backgrnd,self.outdir)
			print(cmd);os.system(cmd)
			self.backgrnd = '%s/%s' % (self.outdir,os.path.basename(self.backgrnd))

		if os.path.splitext(self.phafile)[-1] == '.pha':
			self.basename = '%s_s%dn%d' % (
				os.path.basename(self.phafile).replace('.pha',''),
				self.binminsig, self.binmaxbin)
		elif os.path.splitext(self.phafile)[-1] == '.pi':
			self.basename = '%s_s%dn%d' % (
				os.path.basename(self.phafile).replace('.pi',''),
				self.binminsig, self.binmaxbin)		

		print(self.parerrnum,type(self.parerrnum))

	def get_phafile_property(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.hdu_pha = fits.open(self.phafile)
		self.header = self.hdu_pha['SPECTRUM'].header
		try:
			self.OBSID = self.header['OBS_ID']		
		except:
			self.OBSID = None
		try:
			self.DATEOBS = self.header['DATE-OBS']
		except:			
			self.DATEOBS = None
		try:
			self.DATEEND = self.header['DATE-END']		
		except:
			self.DATEEND = None
		self.EXPOSURE = self.header['EXPOSURE']
		try:
			self.OBJECT = self.header['OBJECT']
		except:
			self.OBJECT = None
		if self.rmffile == None:
			self.rmffile = self.header['RESPFILE']
		if self.arffile == None:
			self.arffile = self.header['ANCRFILE']
		self.TELESCOP = self.header['TELESCOP']
		self.INSTRUME = self.header['INSTRUME']
		try:		
			self.MJDOBS = self.header['MJD-OBS']
		except:
			self.MJDOBS = None
		try:
			self.TSTART = self.header['TSTART']
			self.TSTOP = self.header['TSTOP']	
		except:
			self.TSTART = None
			self.TSTOP = None

		try:
			self.MJD_DATEOBS = float(asttime.Time(self.DATEOBS,format='isot', scale='utc').mjd)
			self.MJD_DATEEND = float(asttime.Time(self.DATEEND,format='isot', scale='utc').mjd)
		except:
			self.MJD_DATEOBS = None
			self.MJD_DATEEND = None

	def show_property(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		#for key, value in self.items(): # python 3 
		#for key, value in vars(self).iteritems(): # python 2 
		for (key, value) in vars(self).items(): # python 2 and 3
			if not key in ["hdu_pha","header"]:
				print(key,value)

	def get_rate_and_error(self,emin,emax):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		tmp_log = '%s/tmp_rate_%sto%skeV.log' % (self.outdir,
			str(emin).replace('.','p'),str(emax).replace('.','p'))
		cmd = 'rm -f %s' % tmp_log
		print(cmd);os.system(cmd)

		cmd  = 'xspec<<EOF\n'
		cmd += 'setplot device tmp/null\n'		
		cmd += 'data 1 %s\n' % self.phafile
		if self.backgrnd is not None:
			cmd += 'back 1 %s\n' % self.backgrnd
		cmd += 'resp 1 %s\n' % self.rmffile
		cmd += 'arf 1 %s\n' % self.arffile
		cmd += 'setplot energy\n'
		cmd += 'ignore 1:**-%0.3f %0.3f-**\n' % (emin,emax)
		cmd += 'log %s\n' % tmp_log
		cmd += 'show rate\n'
		cmd += 'log none\n'
		cmd += 'exit\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)

		rate = grep(tmp_log,"#Net count rate",6)[-1]
		rate_error = grep(tmp_log,"#Net count rate",8)[-1]
	
		#cmd = 'rm -f %s' % tmp_log
		print(cmd);os.system(cmd)	

		outlist = [rate, rate_error]
		return outlist 

	def get_rate_list(self):
		self.ratelist = []
		for i in range(len(self.ratebands)):
			emin = self.ratebands[i][0]
			emax = self.ratebands[i][1]
			rate, rate_error = self.get_rate_and_error(emin,emax)
			self.ratelist.append([rate,rate_error])

			keyword = "rate_%sto%skeV" % (str(emin).replace('.','p'),str(emax).replace('.','p'))
			setattr(self, keyword, rate)
			keyword = "rate_%sto%skeV_err" % (str(emin).replace('.','p'),str(emax).replace('.','p'))
			setattr(self, keyword, rate_error)			

	def bin_spec(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)


		grp_qdp = '%s.qdp' % self.basename
		grp_pco = '%s.pco' % self.basename		
		grp_pha = '%s.pha' % self.basename

		cmd  = 'xspec <<EOF\n'
		cmd += 'setplot device tmp/null\n'
		cmd += '\n'
		cmd += 'data 1 %s\n' % self.phafile
		if self.backgrnd != None:
			cmd += 'back 1 %s\n' % self.backgrnd
		cmd += 'setplot energy\n'
		cmd += 'ignore 1:**-%.2f %.2f-**\n' % (self.fitemin,self.fitemax)
		cmd += 'setplot rebin %d %d\n' % (self.binminsig,self.binmaxbin)
		cmd += 'setplot channel\n'
		cmd += 'iplot data\n'
		cmd += ''
		cmd += 'we %s\n' % self.basename
		cmd += 'exit\n'
		cmd += 'exit\n'
		print(cmd);os.system(cmd)

		cmd  = 'fgrppha.py %s %s' % (self.phafile,grp_qdp)
		print(cmd);os.system(cmd)

		cmd  = 'mv %s %s %s; rm -f %s' % (grp_qdp,grp_pha,self.outdir,grp_pco)
		print(cmd);os.system(cmd)

		self.grp_pha = '%s/%s' % (self.outdir,grp_pha)
		self.phafile = self.grp_pha

	def prepare_read_xcm(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.readxcm = '%s/%s_read.xcm' % (self.outdir,os.path.basename(self.phafile).replace('.pha',''))

		f = open(self.readxcm,'w')
		dump  = 'setplot device tmp/null\n'
		dump += 'data 1 %s\n' % self.phafile	
		if self.backgrnd != None:
			dump += 'back 1 %s\n' % self.backgrnd
		dump += 'resp 1 %s\n' % self.rmffile
		dump += 'arf 1 %s\n' % self.arffile
		dump += 'setplot energy\n'
		dump += 'ignore 1:**-%.2f %.2f-**\n' % (self.fitemin,self.fitemax)
		if self.modelxcm != None:
			dump += '@%s\n' % self.modelxcm
		f.write(dump)
		f.close()

	def make_template_pcofile(self):
		self.fpco_fit = '%s/%s_fit.pco' % (self.outdir,self.basename)
		f = open(self.fpco_fit,'w')
		dump  = 'time off\n'
		dump += 'lwid 5 \n'
		dump += 'lwid 5 on 1..100 \n'	
		dump += 'csize 1.1\n'
		dump += 'lab pos y 2.8\n'
		dump += 'r x %.1f %.1f\n' % (self.plotxmin,self.plotxmax)
		dump += 'r y1 %.1e %.1e\n' % (self.plotymin,self.plotymax)
		dump += 'r y2 %.1f %.1f\n' % (self.ploty2min,self.ploty2max)
		dump += 'col 2 on 2\n'
		dump += 'win 2\n'
		dump += 'LAB  2 COL 2 LIN 0 100 JUS Lef POS %.3f 0 " "\n' % self.fitemin
		f.write(dump)
		f.close()

	def make_template_pcofile_eeuf(self):
		self.fpco_fit_eeuf = '%s/%s_fit_eeuf.pco' % (self.outdir,self.basename)
		f = open(self.fpco_fit_eeuf,'w')
		dump  = 'time off\n'
		dump += 'lwid 5 \n'
		dump += 'lwid 5 on 1..100 \n'	
		dump += 'csize 1.1\n'
		dump += 'lab pos y 2.8\n'
		dump += 'r x %.1f %.1f\n' % (self.plotxmin,self.plotxmax)
		dump += 'r y1 %.1e %.1e\n' % (self.plotymineeuf,self.plotymaxeeuf)
		dump += 'r y2 %.1f %.1f\n' % (self.ploty2min,self.ploty2max)
		dump += 'col 2 on 2\n'
		dump += 'win 2\n'
		dump += 'LAB  2 COL 2 LIN 0 100 JUS Lef POS %.3f 0 " "\n' % self.fitemin
		f.write(dump)
		f.close()

	def fit_spectrum(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.basename_fit = '%s/%s_fit' % (self.outdir,self.basename)
		self.fps_fit = '%s/%s_fit.ps' % (self.outdir,self.basename)
		self.fps_fit_eeuf = '%s/%s_fit_eeuf.ps' % (self.outdir,self.basename)		
		self.fxcm_fit = '%s/%s_fit.xcm' % (self.outdir,self.basename)		
		self.flog_fit = '%s/%s_fit.log' % (self.outdir,self.basename)				

		self.title = '%s %s/%s %s %s (%.1f s)' % (self.OBJECT,
			self.TELESCOP,self.INSTRUME, self.OBSID,
			self.DATEOBS,self.EXPOSURE)
		self.subtitle = '%s' % (os.path.basename(self.phafile))

		cmd  = 'xspec<<EOF\n'
		cmd += '@%s\n' % self.readxcm 
		cmd += 'renorm \n'
		cmd += 'query yes\n'	
		cmd += 'fit\n'		
		cmd += 'log %s\n' % self.flog_fit 
		cmd += 'save all %s\n' % self.fxcm_fit						
		cmd += 'show rate\n'
		cmd += 'show pa\n'		
		cmd += 'show fit\n'	
		cmd += 'log none\n'
		cmd += 'iplot ld del\n'
		if os.path.exists(self.fpco_fit):
			cmd += '@%s\n' % self.fpco_fit
		cmd += 'win 1 \n'
		cmd += 'la t %s\n' % self.title
		cmd += 'la f %s\n' % self.subtitle
		cmd += 'hard %s/cps\n' % self.fps_fit  
		#cmd += 'we %s\n' % self.basename_fit
		cmd += 'we temp\n' 
		cmd += 'y\n'
		cmd += 'exit\n'
		cmd += 'iplot eeuf del\n'
		if os.path.exists(self.fpco_fit_eeuf):
			cmd += '@%s\n' % self.fpco_fit_eeuf	
		cmd += 'win 1 \n'
		cmd += 'la t %s\n' % self.title
		cmd += 'la f %s\n' % self.subtitle				
		cmd += 'hard %s/cps\n' % self.fps_fit_eeuf 
		#cmd += 'we %s\n' % self.basename_fit
		cmd += 'we temp_eeuf\n' 
		cmd += 'y\n'
		cmd += 'exit\n'		
		cmd += 'exit\n'
		cmd += 'EOF\n'			
		print(cmd);os.system(cmd)

		#cmd  =  "grep -l 'temp.pco' temp.qdp | xargs sed -i.bak -e 's/temp.pco/%s.pco/g';" % os.path.basename(self.basename_fit)
		cmd  =  "grep -l 'temp.pco' temp.qdp | xargs sed -e 's/temp.pco/%s.pco/g';" % os.path.basename(self.basename_fit)
		cmd += 'mv temp.qdp %s.qdp;' % self.basename_fit
		cmd += 'mv temp.pco %s.pco;' % self.basename_fit		
		print(cmd);os.system(cmd)

		cmd = 'ps2pdf.py %s' % self.fps_fit
		print(cmd);os.system(cmd)

		cmd  =  "grep -l 'temp_eeuf.pco' temp_eeuf.qdp | xargs sed -e 's/temp_eeuf.pco/%s_eeuf.pco/g';" % os.path.basename(self.basename_fit)
		cmd += 'mv temp_eeuf.qdp %s_eeuf.qdp;' % self.basename_fit
		cmd += 'mv temp_eeuf.pco %s_eeuf.pco;' % self.basename_fit		
		print(cmd);os.system(cmd)

		cmd = 'ps2pdf.py %s' % self.fps_fit_eeuf
		print(cmd);os.system(cmd)		

	def get_flux(self,fitxcm,emin,emax):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		tmp_log = '%s/tmp_flux_%sto%skeV.log' % (self.outdir,
			str(emin).replace('.','p'),str(emax).replace('.','p'))
		cmd = 'rm -f %s' % tmp_log
		print(cmd);os.system(cmd)

		cmd  = 'xspec<<EOF\n'
		cmd += '@%s\n' % self.fxcm_fit
		cmd += 'query yes\n'		
		cmd += 'fit\n'
		cmd += 'log %s\n' % tmp_log
		cmd += 'flux %.3f %.3f err 300,68.3\n' % (emin,emax)
		cmd += 'log none\n'
		cmd += 'exit\n'
		print(cmd);os.system(cmd)

		for line in open(tmp_log):
			cols = line.split()
			if "# Model Flux" in line:
				flux = float(cols[5].replace('(',''))
			if "Error range" in line:
				flux_min = float(cols[6].replace('(',''))
				flux_max = float(cols[8].replace(')',''))				
		flux_error_min = flux - flux_min
		flux_error_max = flux_max - flux
		outlist = [flux,flux_error_min,flux_error_max]
		return outlist 

	def get_flux_list(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.fluxlist = []
		for i in range(len(self.fluxbands)):
			emin = self.fluxbands[i][0]
			emax = self.fluxbands[i][1]
			flux, flux_error_min, flux_error_max = self.get_flux(self.fxcm_fit,emin,emax)
			self.fluxlist.append([flux, flux_error_min, flux_error_max])

			keyword = "flux_%sto%skeV" % (str(emin).replace('.','p'),str(emax).replace('.','p'))
			setattr(self, keyword, flux)
			keyword = "flux_%sto%skeV_err_min" % (str(emin).replace('.','p'),str(emax).replace('.','p'))
			setattr(self, keyword, flux_error_min)	
			keyword = "flux_%sto%skeV_err_max" % (str(emin).replace('.','p'),str(emax).replace('.','p'))
			setattr(self, keyword, flux_error_max)				

	def get_parerror(self,fitxcm,parnum):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		tmp_log = '%s/tmp_par%d.log' % (self.outdir,parnum)
		cmd = 'rm -f %s' % tmp_log
		print(cmd);os.system(cmd)

		cmd  = 'xspec<<EOF\n'
		cmd += '@%s\n' % self.fxcm_fit
		cmd += 'query yes\n'
		cmd += 'fit\n'
		cmd += 'log %s\n' % tmp_log
		cmd += 'error 1.0 %d\n' % parnum
		cmd += 'log none\n'
		cmd += 'exit\n'
		print(cmd);os.system(cmd)

		print(tmp_log)

		flag_show_pa = False
		for line in open(self.flog_fit):
			cols = line.split()
			if flag_show_pa and len(cols) > 1:
				if str(cols[1]) == str(parnum):
					value = float(cols[-3])
					break
			if "#Parameters defined:" in line:
				flag_show_pa = True
		if not flag_show_pa:
			value = np.nan				


		flag_err = False
		for line in open(tmp_log):
			cols = line.split()
			if flag_err and len(cols) == 5:
				if cols[1] == str(parnum) and cols[4][0] == '(' and cols[4][-1] == ')':
					print(line)
					tmp1,tmp2 = cols[4].split(',')
					err_min = float(tmp1.replace('(',''))
					err_max = float(tmp2.replace(')',''))
					break
			if "# Parameter   Confidence Range (1)" in line:
				flag_err = True
		if not flag_err:
			err_min = np.nan
			err_max = np.nan
	
		outlist = [value,err_min,err_max]
		return outlist 		

	def get_parerror_list(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.parerrorlist = []
		for parnum in self.parerrnum:
			value, err_min, err_max = self.get_parerror(self.fxcm_fit,parnum)
			self.parerrorlist.append([value, err_min, err_max])

			keyword = "par_%d" % parnum
			setattr(self, keyword, value)
			keyword = "par_%d_err_min" % parnum
			setattr(self, keyword, err_min)	
			keyword = "par_%d_err_max" % parnum
			setattr(self, keyword, err_max)				

	def dump_yamlfile(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		param = {}
		#for key, value in vars(self).iteritems(): # python 2 
		for (key, value) in vars(self).items(): # python 2 and 3
			if not key in ["hdu_pha","header","fitlog"]:
				param[key] = value 
				print(key,value)
	
		#self.fyaml_dump = '%s/%s.yaml' % (self.outdir,self.basename)
		self.fyaml_dump = '%s/%s_fit.yaml' % (self.outdir,self.basename)
		with open(self.fyaml_dump, "w") as wf:
		    yaml.dump(param, wf,default_flow_style=False)		

	def analyze_fitlog(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.fitlog = XspecFitLog(self.flog_fit)		
		self.chisquare = self.fitlog.get_chisquare()
		self.dof = self.fitlog.get_dof()
		self.reduced_chisquare = self.fitlog.get_reduced_chisquare()		
		self.probability = self.fitlog.get_probability()

		print(self.chisquare)
		print(self.reduced_chisquare)
		print(self.dof)
		print(self.probability)

	def set_flag_complete(self):
		self.flag_fit_complete = True

	def run(self):
		self.get_phafile_property()
		self.get_rate_list()
		self.bin_spec()
		self.show_property()
		self.prepare_read_xcm()
		self.make_template_pcofile()
		self.make_template_pcofile_eeuf()		
		self.fit_spectrum()
		self.analyze_fitlog()
		self.get_flux_list()
		self.get_parerror_list()
		self.set_flag_complete()
		self.dump_yamlfile()

class XspecFitLog():
	def __init__(self, logfile):
		self.logfile = logfile
		if not os.path.exists(self.logfile):
			print("file does not exits: %s" % self.logfile)
			quit()

	def get_chisquare(self):
		#self.chisquare = grep(self.logfile,"#Test statistic : Chi-Squared =",5)[-1]
		self.chisquare = grep(self.logfile,"#Test statistic :",4)[-1]
		return self.chisquare

	def get_reduced_chisquare(self):
		#self.reduced_chisquare = grep(self.logfile,"# Reduced chi-squared =",4)[-1]
		self.reduced_chisquare = self.chisquare / float(self.dof)
		return self.reduced_chisquare

	def get_dof(self):
		#self.dof = grep(self.logfile,"# Reduced chi-squared =",6,dtype=int)[-1]
		self.dof = grep(self.logfile,"# Null hypothesis probability of",7,dtype=int)[-1]
		return self.dof		

	def get_probability(self):
		#self.probability = grep(self.logfile,"# Null hypothesis probability =",5)[-1]		
		self.probability = grep(self.logfile,"# Null hypothesis probability of",5,dtype=float)[-1]				
		return self.probability		

class CSVtoXSPEC():
	def __init__(self,filelist,yamlfile):
		self.filelist = filelist
		self.yamlfile = yamlfile

		if not os.path.exists(self.filelist):
			sys.stderr.write('file %s does not exist.' % self.filelist)
			exit()
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('file %s does not exist.' % self.yamlfile)
			exit()			

		self.param = yaml.load(open(self.yamlfile),Loader=yaml.FullLoader)

	def make_csv2xspec(self,outcsvfile):
		df = pd.read_csv(self.filelist,names=('data_id','phafile','backgrnd','rmffile','arffile'))
		df["modelxcm"] = [self.param["modelxcm"] for i in range(len(df))]
		df["binminsig"] = [self.param["binminsig"] for i in range(len(df))]
		df["binmaxbin"] = [self.param["binmaxbin"] for i in range(len(df))]		
		df["fitemin"] = [self.param["fitemin"] for i in range(len(df))]				
		df["fitemax"] = [self.param["fitemax"] for i in range(len(df))]						
		df["plotxmin"] = [self.param["plotxmin"] for i in range(len(df))]								
		df["plotxmax"] = [self.param["plotxmax"] for i in range(len(df))]										
		df["plotymin"] = [self.param["plotymin"] for i in range(len(df))]												
		df["plotymax"] = [self.param["plotymax"] for i in range(len(df))]
		df["plotymineeuf"] = [self.param["plotymineeuf"] for i in range(len(df))]												
		df["plotymaxeeuf"] = [self.param["plotymaxeeuf"] for i in range(len(df))]																														
		df["ploty2min"] = [self.param["ploty2min"] for i in range(len(df))]																
		df["ploty2max"] = [self.param["ploty2max"] for i in range(len(df))]																		
		#df["flagrun"] = [1 for i in range(len(df))]																				
		print(df)

		outdir = os.path.dirname(outcsvfile)
		if not os.path.exists(outdir):
			cmd = 'mkdir -p %s' % outdir
			print(cmd);os.system(cmd)

		df.to_csv(outcsvfile)


class MonitoringManager():
	def __init__(self,csvfile,yamlfile):
		self.csvfile = csvfile 
		self.yamlfile = yamlfile

		if not os.path.exists(self.csvfile):
			sys.stderr.write('error: file %s does not exist.\n' % self.csvfile)
			exit()			
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('error: file %s does not exist.\n' % self.yamlfile)
			exit()						
		print("csvfile: %s" % self.csvfile)
		print("yamlfile: %s" % self.yamlfile)

		self.df = pd.read_csv(self.csvfile)
		self.df = self.df.drop(columns='Unnamed: 0')
		self.param = yaml.load(open(self.yamlfile),Loader=yaml.FullLoader)

		self.outcsvfile = '%s/%s' % (
			self.param['outdir'],
			os.path.basename(self.csvfile).replace('.csv','_fit.csv'))
		print("out csvfile: %s" % self.outcsvfile)

	def run_fit(self):
		for index, dataset in self.df.iterrows():
			outdir = '%s/%s' % (self.param['outdir'],dataset['data_id'])
			xspec_pha = XspecPha(
				phafile=dataset['phafile'],
				outdir=outdir,
				backgrnd=dataset['backgrnd'],
				rmffile=dataset['rmffile'],
				arffile=dataset['arffile'],
				modelxcm=dataset['modelxcm'],
				binminsig=dataset['binminsig'],binmaxbin=dataset['binmaxbin'],
				fitemin=dataset['fitemin'],fitemax=dataset['fitemax'],
				plotxmin=dataset['plotxmin'],plotxmax=dataset['plotxmax'],
				plotymin=dataset['plotymin'],plotymax=dataset['plotymax'],
				plotymineeuf=dataset['plotymineeuf'],plotymaxeeuf=dataset['plotymaxeeuf'],				
				ploty2min=dataset['ploty2min'],ploty2max=dataset['ploty2max'],
				ratebands=self.param['ratebands'],
				fluxbands=self.param['fluxbands'],
				parerrnum=self.param['parerrnum'])
			xspec_pha.run()			
			#try:				
			#	xspec_pha.run()
			#except:
			#	print("can not fit ... skip %s" % outdir)

	def collect_fitresults(self):	
		add_column_names = [
			"grp_pha",
			"OBSID",
			"DATEOBS",
			"DATEEND",
			"EXPOSURE",
			"MJD_DATEOBS",
			"MJD_DATEEND",
			"OBJECT",
			"TELESCOP",
			"INSTRUME",
			"MJDOBS",
			"TSTART",
			"TSTOP",
			"fxcm_fit",
			"title",
			"chisquare",
			"reduced_chisquare",
			"dof",
			"probability"
			]

		for ebands in self.param['ratebands']:
			add_column_names.append("rate_%sto%skeV" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))
			add_column_names.append("rate_%sto%skeV_err" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))			

		for ebands in self.param['fluxbands']:
			add_column_names.append("flux_%sto%skeV" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))
			add_column_names.append("flux_%sto%skeV_err_min" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))			
			add_column_names.append("flux_%sto%skeV_err_max" % (str(ebands[0]).replace('.','p'),str(ebands[1]).replace('.','p')))						

		for parnum in self.param['parerrnum']:
			add_column_names.append("par_%d" % parnum)
			add_column_names.append("par_%d_err_min" % parnum)
			add_column_names.append("par_%d_err_max" % parnum)

		add_dictionary = {}
		for keyword in add_column_names:
			add_dictionary[keyword] = []

		for index, dataset in self.df.iterrows():
			basename = os.path.splitext(os.path.basename(dataset['phafile']))[0]
			fitresult_yamlfile = glob.glob('%s/%s/%s*.yaml' % (self.param['outdir'],dataset['data_id'],basename))[0]
			fitresult = yaml.load(open(fitresult_yamlfile),Loader=yaml.FullLoader)
			print("\n")
			print(fitresult)

			for keyword in add_column_names:
				print(keyword)
				add_dictionary[keyword].append(fitresult[keyword])
		
		self.df_add = pd.DataFrame(add_dictionary,columns=(add_column_names))
		self.df_new = self.df.join([self.df_add])
		self.df_new.to_csv(self.outcsvfile)

