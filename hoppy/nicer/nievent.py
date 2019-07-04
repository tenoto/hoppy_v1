#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys 
import yaml

class NicerEventFits():
	def __init__(self,inputfits):
		self.inputfits = inputfits
		if not os.path.exists(self.inputfits):
			sys.stderr.write('file %s does not exist.' % self.inputfits)
			exit()			

	def barycentric_correction(self,ra,dec,orbfile,
		outfits=None,refframe="ICRS",ephem="JPLEPH.430"):
		print("--- barycentric_correction ---")	

		if not os.path.exists(orbfile):
			sys.stderr.write('file %s does not exist.' % orbfile)
			exit()		
		if outfits != None:
			self.baryfits = outfits 
			self.baryfits_log = self.baryfits.replace('.evt','.log')
		else:
			self.baryfits = os.path.basename(self.inputfits.replace('.gz','').replace('.evt','_bary.evt'))
			self.baryfits_log = self.baryfits.replace('.evt','.log')

		if not os.path.exists(os.path.dirname(self.baryfits)):
			cmd = 'mkdir -p %s' % os.path.dirname(self.baryfits)
			print(cmd);os.system(cmd)

		print("output file: %s " % self.baryfits)
		print("orbit file: %s" % orbfile)
		print("RA (J2000): %.6f" % ra)
		print("DEC (J2000): %.6f" % dec)		
		print("refframe: %s" % refframe)				
		print("ephem: %s" % ephem)						

		cmd = 'rm -f %s' % self.baryfits
		print(cmd);os.system(cmd)

		cmd  = "barycorr infile=%s " % self.inputfits
		cmd += "outfile=%s " % self.baryfits
		cmd += "ra=%.6f dec=%.6f " % (ra,dec)
		cmd += "orbitfiles=%s " % orbfile
		cmd += "refframe=%s ephem=%s " % (refframe,ephem)
		cmd += ">& %s" % self.baryfits_log
		print(cmd);os.system(cmd)

		cmd  = 'fcalc clobber=yes infile=%s[\'EVENTS\'] ' % self.baryfits
		cmd += 'outfile=%s clname="BARY_TIME" ' % self.baryfits
		cmd += 'expr="#MJDREFI + #MJDREFF + TIME/86400.0" '
		print(cmd);os.system(cmd)

		cmd  = 'fcalc clobber=yes infile=%s[\'GTI\'] ' % self.baryfits
		cmd += 'outfile=%s clname="BARY_START" ' % self.baryfits
		cmd += 'expr="#MJDREFI + #MJDREFF + START/86400.0" '
		print(cmd);os.system(cmd)

		cmd  = 'fcalc clobber=yes infile=%s[\'GTI\'] ' % self.baryfits
		cmd += 'outfile=%s clname="BARY_STOP" ' % self.baryfits
		cmd += 'expr="#MJDREFI + #MJDREFF + STOP/86400.0" '
		print(cmd);os.system(cmd)		

	"""
	def run_photonphase(self,fyaml,orbfile):
		self.fyaml = fyaml
		print("fyaml : %s" % fyaml)
		if not os.path.exists(self.fyaml):
			sys.stderr.write('file %s does not exist.' % self.fyaml)
			exit()
		self.param = yaml.load(open(self.fyaml))

		print("orbfile : %s" % orbfile)			
		if not os.path.exists(orbfile):
			sys.stderr.write('file %s does not exist.' % orbfile)
			exit()					

		self.barycentric_correction(
			self.param['RA'],self.param['DEC'],orbfile,outfits=None,
			refframe=self.param['REFFRAME'],ephem=self.param['EPHEM'])

	def addphase(self):
		self.outphaseevt = self.outbaryevt.replace('.evt','_phase.evt')
		cmd = 'rm -f %s' % self.outphaseevt
		print(cmd);os.system(cmd)
				
		cmd  = "faddphase_nu.py "
		cmd += "-i %s " % self.outbaryevt
		cmd += "-o %s " % self.outphaseevt 
		cmd += "--nu %.15f " % self.param['nu']
		cmd += "--nudot %.6e " % self.param['nudot']
		cmd += "-e %.6f " % self.param['epoch']
		print(cmd);os.system(cmd)

		cmd = 'rm -f %s' % self.outbaryevt
		print(cmd);os.system(cmd)

	def plot_profile(self):
		self.profile = self.outphaseevt.replace('.evt','_pls.fits')
		cmd = 'rm -f %s' % self.profile
		print(cmd);os.system(cmd)

		cmd = 'rm -f %s' % self.profile.replace('.fits','.pdf')
		print(cmd);os.system(cmd)		

		cmd  = 'fplot_pulseprofile.py '
 		cmd += '-i %s ' % self.outphaseevt 
 		cmd += '-o %s ' % self.profile
 		cmd += '-n %d ' % self.param['phase_bin']
		cmd += '-c PHASE'
		print(cmd);os.system(cmd)

	def addabsphase(self):
		self.outphaseevt2 = self.outphaseevt.replace('.evt','2.evt')

		cmd  = "faddphase_nu.py "
		cmd += "-i %s " % self.outphaseevt 
		cmd += "-o %s " % self.outphaseevt2 
		cmd += "--nu %.15f " % self.param['nu']
		cmd += "--nudot %.6e " % self.param['nudot']
		cmd += "-e %.6f " % self.param['epoch']
		cmd += "--absphase "
		print(cmd);os.system(cmd)

		cmd = 'rm -f %s' % self.outphaseevt
		print(cmd);os.system(cmd)		
	"""		