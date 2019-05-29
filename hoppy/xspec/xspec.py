# -*- coding: utf-8 -*-

import os 
import sys 
import yaml 
import argparse 
import numpy as np
import astropy.io.fits as pyfits 
import astropy.time as asttime 

def string_to_list(string):
	mainlist = []
	for col in string.split(','):
		sublist = []
		for i in col.split('-'):
			sublist.append(float(i))
		mainlist.append(sublist)
	return mainlist

class XspecPha():
	def __init__(self,phafile,outdir='out',
		backgrnd=None,rmffile=None,arffile=None,modelxcm=None,
		binminsig=5,binmaxbin=50,fitemin=0.4,fitemax=10.0,
		plotxmin=0.4,plotxmax=10.0,plotymin=1e-4,plotymax=1e+4,
		plotymineeuf=1e-4,plotymaxeeuf=10,ploty2min=-10.0,ploty2max=10.0,
		ratebands=[[0.4,6.0],[1.0,10.0]],
		fluxbands=[[0.4,6.0],[1.0,10.0]],
		parerrnum=[1,2,5]):
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
		self.ratebands = ratebands
		self.fluxbands = fluxbands
		self.parerrnum = parerrnum

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
		cmd = 'cp %s %s' % (self.backgrnd,self.outdir)
		print(cmd);os.system(cmd)
		self.backgrnd = '%s/%s' % (self.outdir,os.path.basename(self.backgrnd))

	def get_phafile_property(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.hdu_pha = pyfits.open(self.phafile)
		self.header = self.hdu_pha['SPECTRUM'].header
		self.OBSID = self.header['OBS_ID']		
		self.DATEOBS = self.header['DATE-OBS']
		self.DATEEND = self.header['DATE-END']		
		self.EXPOSURE = self.header['EXPOSURE']
		self.OBJECT = self.header['OBJECT']
		if self.rmffile == None:
			self.rmffile = self.header['RESPFILE']
		if self.arffile == None:
			self.arffile = self.header['ANCRFILE']
		self.TELESCOP = self.header['TELESCOP']
		self.INSTRUME = self.header['INSTRUME']		
		self.MJDOBS = self.header['MJD-OBS']
		self.TSTART = self.header['TSTART']
		self.TSTOP = self.header['TSTOP']	

		self.MJD_DATEOBS = asttime.Time(self.DATEOBS,format='isot', scale='utc').mjd	
		self.MJD_DATEEND = asttime.Time(self.DATEEND,format='isot', scale='utc').mjd

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
			self.ratelist.append(self.get_rate_and_error(emin,emax))

	def bin_spec(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.basename = '%s_s%dn%d' % (
			os.path.basename(self.phafile).replace('.pha',''),
			self.binminsig, self.binmaxbin)
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
			self.fluxlist.append(self.get_flux(self.fxcm_fit,emin,emax))

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

		flag_show_pa = False
		for line in open(self.flog_fit):
			cols = line.split()
			if flag_show_pa and len(cols) > 1:
				if cols[1] == str(parnum):
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
			self.parerrorlist.append(self.get_parerror(self.fxcm_fit,parnum))

	def dump_yamlfile(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		param = {}
		#for key, value in vars(self).iteritems(): # python 2 
		for (key, value) in vars(self).items(): # python 2 and 3
			if not key in ["hdu_pha","header","fitlog"]:
				param[key] = value 
				print(key,value)
	
		self.fyaml_dump = '%s/%s.yaml' % (self.outdir,self.basename)
		with open(self.fyaml_dump, "w") as wf:
		    yaml.dump(param, wf,default_flow_style=False)		

	def analyze_fitlog(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.fitlog = XspecFitLog(self.flog_fit)		
		self.chisquare = self.fitlog.get_chisquare()
		self.reduced_chisquare = self.fitlog.get_reduced_chisquare()
		self.dof = self.fitlog.get_dof()
		self.probability = self.fitlog.get_probability()

		print(self.chisquare)
		print(self.reduced_chisquare)
		print(self.dof)
		print(self.probability)

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
		self.dump_yamlfile()

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

class XspecFitLog():
	def __init__(self, logfile):
		self.logfile = logfile
		if not os.path.exists(self.logfile):
			print("file does not exits: %s" % self.logfile)
			quit()

	def get_chisquare(self):
		self.chisquare = grep(self.logfile,"#Test statistic : Chi-Squared =",5)[-1]
		return self.chisquare

	def get_reduced_chisquare(self):
		self.reduced_chisquare = grep(self.logfile,"# Reduced chi-squared =",4)[-1]
		return self.reduced_chisquare

	def get_dof(self):
		self.dof = grep(self.logfile,"# Reduced chi-squared =",6,dtype=int)[-1]
		return self.dof		

	def get_probability(self):
		self.probability = grep(self.logfile,"# Null hypothesis probability =",5)[-1]
		return self.probability		



	#def get_fitrange_rate(logfile):
	#	for line in open(logfile):
	#		cols = line.split()
	#		if "#Net count rate" in line:
	#			rate = float(cols[6])
	#			rate_err = float(cols[8])
	#			break
	#	return rate, rate_err			

if __name__=="__main__":

	sys.stdout.write('\n... run a single xspec fitting ...\n')

	parser = argparse.ArgumentParser(
		prog=__file__,
		usage='python %s phafile [-b backgrnd] [-r rmffile] [-a arffile] [-m modelxcm] [-s binminsig] [-n binmaxbin]' % os.path.basename(__file__),
		description='Automatic Xspec fitting.',
		epilog='',
		add_help=True,
		)
	parser.add_argument(
		'phafile',metavar='phafile',type=str,
		help='source pha file for fitting.') 
	parser.add_argument(
		'-o','--outdir',metavar='outdir',type=str,default='out',
		help='output directory.') 		
	parser.add_argument(
		'-b','--backgrnd',metavar='backgrnd',type=str,default=None,
		help='background pha file for fitting.') 		
	parser.add_argument(
		'-r','--rmffile',metavar='rmffile',type=str,default=None,
		help='rmffile for fitting.') 
	parser.add_argument(
		'-a','--arffile',metavar='arffile',type=str,default=None,
		help='arffile for fitting.') 
	parser.add_argument(
		'-m','--modelxcm',metavar='modelxcm',type=str,default=None,
		help='modelxcm for fitting.') 				
	parser.add_argument(
		'-s','--binminsig',metavar='binminsig',type=int,default=5,
		help='binning minimum significance.') 					
	parser.add_argument(
		'-n','--binmaxbin',metavar='binmaxbin',type=int,default=50,
		help='binning max bins.') 
	parser.add_argument(
		'--fitemin',metavar='fitemin',type=float,default=0.4,
		help='fitting energy min (keV).') 								
	parser.add_argument(
		'--fitemax',metavar='fitemax',type=float,default=10.0,
		help='fitting energy max (keV).') 	
	parser.add_argument(
		'--ratebands',metavar='ratebands',type=str,default=[[0.4,6.0],[1.0,10.0]],
		help='rate energy bands (list) example:0.8-6.0,2.0-10.0 .') 											
	parser.add_argument(
		'--fluxbands',metavar='fluxbands',type=str,default=[[0.4,6.0],[1.0,10.0]],
		help='flux energy bands (list) example:0.8-6.0,2.0-10.0 .') 	
	parser.add_argument(
		'--parerrnum',metavar='parerrnum',type=str,default=[1,2,5],
		help='parameter error number list.') 																
	args = parser.parse_args()	
	print(args)

	if type(args.ratebands) == list:
		ratebands = args.ratebands
	elif type(args.ratebands) == str:
		ratebands = string_to_list(args.ratebands)

	if type(args.fluxbands) == list:
		fluxbands = args.fluxbands
	elif type(args.fluxbands) == str:
		fluxbands = string_to_list(args.fluxbands)		

	if type(args.parerrnum) == list:
		parerrnum = args.parerrnum
	elif type(args.fluxbands) == str:
		parerrnum = string_to_list(args.parerrnum)	

	xspec_pha = XspecPha(args.phafile,
		outdir=args.outdir,
		backgrnd=args.backgrnd,rmffile=args.rmffile,arffile=args.arffile,modelxcm=args.modelxcm,
		binminsig=args.binminsig,binmaxbin=args.binmaxbin,fitemin=args.fitemin,fitemax=args.fitemax,
		ratebands=ratebands,fluxbands=fluxbands,parerrnum=parerrnum)
	xspec_pha.run()


