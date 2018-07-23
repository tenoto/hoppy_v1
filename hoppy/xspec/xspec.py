# -*- coding: utf-8 -*-

import os 
import sys 
import argparse 
import astropy.io.fits as pyfits 

class XspecPha():
	def __init__(self,phafile,outdir='out',
		backgrnd=None,rmffile=None,arffile=None,modelxcm=None,
		binminsig=5,binmaxbin=50,fitemin=0.4,fitemax=10.0,
		plotxmin=0.4,plotxmax=10.0,plotymin=1e-4,plotymax=1e+4,
		ploty2min=-10.0,ploty2max=10.0):
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
		self.ploty2min = ploty2min
		self.ploty2max = ploty2max

		if not os.path.exists(self.phafile):
			sys.stderr.write('phafile %s does not exist.' % self.phafile)
			exit()
		if not os.path.exists(self.outdir):
			cmd = 'mkdir -p %s' % self.outdir
			print(cmd);os.system(cmd)

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

	def show_property(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		#for key, value in self.items(): # python 3 
		for key, value in vars(self).iteritems(): # python 2 
			if not key in ["hdu_pha","header"]:
				print(key,value)

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

		self.phafile = '%s/%s' % (self.outdir,grp_pha)

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
		dump += 'r y %.1e %.1e\n' % (self.plotymin,self.plotymax)
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
		self.fxcm_fit = '%s/%s_fit.xcm' % (self.outdir,self.basename)		
		self.flog_fit = '%s/%s_fit.log' % (self.outdir,self.basename)				

		self.title = '%s %s/%s %s' % (self.OBJECT,
			self.TELESCOP,self.INSTRUME, self.OBSID)
		self.subtitle = '%s' % (self.DATEOBS)

		cmd  = 'xspec<<EOF\n'
		cmd += '@%s\n' % self.readxcm 
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
		cmd += 'we %s\n' % self.basename_fit
		cmd += 'y\n'
		cmd += 'exit\n'
		cmd += 'exit\n'
		cmd += 'EOF\n'			
		print(cmd);os.system(cmd)

		cmd = 'ps2pdf.py %s' % self.fps_fit
		print(cmd);os.system(cmd)

	def run(self):
		self.get_phafile_property()
		self.bin_spec()
		self.show_property()
		self.prepare_read_xcm()
		self.make_template_pcofile()
		self.fit_spectrum()

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
	args = parser.parse_args()	
	print(args)

	xspec_pha = XspecPha(args.phafile,
		outdir=args.outdir,
		backgrnd=args.backgrnd,rmffile=args.rmffile,arffile=args.arffile,modelxcm=args.modelxcm,
		binminsig=args.binminsig,binmaxbin=args.binmaxbin,fitemin=args.fitemin,fitemax=args.fitemax)
	xspec_pha.run()


