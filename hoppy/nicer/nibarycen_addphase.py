#!/usr/bin/env python

import os 
import sys 
import yaml 
import argparse
#import hoppy.nicer.nipipeline as nip

class NiPulsar():
	def __init__(self,fyaml):
		if not os.path.exists(fyaml):
			sys.stderr.write('file %s does not exist.' % fyaml)
			exit()
		self.fyaml = fyaml
		print("input yaml file: %s" % self.fyaml)

		self.param = yaml.load(open(self.fyaml))
		print(self.param)

	def barycentric_correction(self):
		if not os.path.exists(self.param['inputevt']):
			sys.stderr.write('file %s does not exist.' % self.param['inputevt'])
			exit()		
		if not os.path.exists(self.param['orbitfile']):
			sys.stderr.write('file %s does not exist.' % self.param['orbitfile'])
			exit()		

		self.outbaryevt = '%s_bary.evt' % os.path.basename(self.param['inputevt']).replace('.evt.gz','').replace('.evt','')
		cmd = 'rm -f %s' % self.outbaryevt
		print(cmd);os.system(cmd)

		cmd  = "barycorr infile=%s " % self.param['inputevt']
		cmd += "outfile=%s " % self.outbaryevt
		cmd += "ra=%.6f dec=%.6f " % (self.param['pulsar_ra'],self.param['pulsar_dec'])
		cmd += "orbitfiles=%s " % self.param['orbitfile']
		cmd += "refframe=ICRS ephem=JPLEPH.430 "
		print(cmd);os.system(cmd)

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

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='nibarycen_addphase.py',
		usage='nibarycen_addphase.py param.yaml',
		description='NICER baryenctric correction and add pulsar phase from a pulsar ephemeris',
		epilog='',
		add_help=True)
	parser.add_argument(
		'fyaml',metavar='fyaml',type=str,        
		help='yaml file with parameters')
	args = parser.parse_args()	

	nipulsar = NiPulsar(args.fyaml)
	nipulsar.barycentric_correction()
	nipulsar.addphase()
	nipulsar.plot_profile()
	nipulsar.addabsphase()
