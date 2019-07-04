#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script add "PULSE_PHASE" column calculated from 
# nu, nudot, nuddot, and epoch. 
# 
# HISTORY
# 2015-02-07: The original perl script aeplsphase.pl is written by Teruaki Enoto.
# 2017-08-26: Translate to python script (Teru Enoto, ver 1.00)
# 2018-09-14: Modified to compatible to photonphase (PINT) (Teru Enoto, ver 2.00)

__author__  = 'Teru Enoto'
__date__    = '2018 September 14'
__version__ = '2.00'

import os 
import sys 
import argparse
from datetime import datetime 
#import hoppy.nicer.nievent as nievt

class EventFits():
	def __init__(self,infits):
		sys.stdout.write('Object of EventFits is created.\n')
		self.infits = infits
		sys.stdout.write('input fits file: %s\n' % self.infits)
		if not os.path.exists(self.infits):
			sys.stderr.write('input file does not exists: %s\n' % self.infits)
			quit()

	def faddphase_nu(self,epoch,
		nu,nudot=0.0,nu2dot=0.0,nu3dot=0.0,nu4dot=0.0,
		outfits=None,offset=0.0,flag_mjd=False):
		sys.stdout.write("--%s--\n" % sys._getframe().f_code.co_name)

		if outfits != None:
			outfits = outfits 			
		else:
			outfits = self.infits.replace('.evt','_phase.evt').replace('.fits','_phase.fits')
		if os.path.exists(outfits):
			sys.stderr.write('error: file %s has already existed.\n' % outfits)
			quit()
		outfits_log = outfits.replace('.evt','.log').replace('.fits','.log')

		if flag_mjd:
			operation1 = "(%.12e*(BARY_TIME-(%.15f))*86400. + (%.7e)*((BARY_TIME-(%.15f))*86400.)**2/2.0 + (%.7e)*((BARY_TIME-(%.15f))*86400.)**3/6.0 + (%.7e)*((BARY_TIME-(%.15f))*86400.)**4/24.0 + (%.7e)*((BARY_TIME-(%.15f))*86400.)**5/120.0)" % (
				nu, epoch, nudot, epoch, nu2dot, epoch, nu3dot, epoch, nu4dot, epoch) 
		else:
			operation1 = "(%.12e*(TIME-(%.7f)) + (%.7e)*(TIME-(%.7f))**2/2.0 + (%.7e)*(TIME-(%.7f))**3/6.0 + (%.7e)*(TIME-(%.7f))**4/24.0 + (%.7e)*(TIME-(%.7f))**5/120.0)" % (
				nu, epoch, nudot, epoch, nu2dot, epoch, nu3dot, epoch, nu4dot, epoch) 
		operation2 = 'floor( PHASE )'		
		operation3 = 'PHASE - PULSE_NUMBER'	
		operation4 = 'PULSE_PHASE < 0.5? PULSE_NUMBER: PULSE_NUMBER+1'
		#operation2 = '(PHASE %% 1.0) + %.7f' % (offset)
		#operation3 = ' PULSE_PHASE<0 ? PULSE_PHASE+1 : PULSE_PHASE '
		#operation4 = 'floor( PHASE )'

		history_dump = """
HISTORY -----------------------------------------------------
HISTORY faddphase_nu version %s at %s 
HISTORY -----------------------------------------------------
HISTORY inputfits ='%s'
HISTORY outputfits='%s'
HISTORY nu = %.12f / Frequency (Hz) used at PULSE_PHAES calcuration
HISTORY nudot  = %.7e   / Frequency derivative (Hz/sec) 
HISTORY nu2dot = %.7e   / Frequency 2nd derivative (Hz/sec2) 
HISTORY nu3dot = %.7e   / Frequency 3rd derivative (Hz/sec3) 
HISTORY nu4dot = %.7e   / Frequency 4th derivative (Hz/sec4) 
HISTORY epoch  = %.7f   / epoch (s) 
HISTORY offset = %.7e   / Optional phase offset 
HISTORY operation1: %s 
HISTORY operation2: %s 
HISTORY operation3: %s 
""" % (__version__, datetime.now().strftime('%Y-%m-%dT%H:%M:%S'),
		self.infits, outfits,
		nu, nudot, nu2dot, nu3dot, nu4dot, epoch, offset, 
		operation1,operation2,operation3)
		print(history_dump)

		cmd  = 'fcalc infile=%s+1 ' % self.infits
		cmd += 'outfile=%s ' % outfits
		cmd += 'clname=\"PHASE\" expr=\"%s\" rowrange=\"-\"' % operation1
		print(cmd); os.system(cmd)

		cmd  = 'fcalc clobber=yes infile=%s+1 ' % outfits 
		cmd += 'outfile=%s ' % outfits
		cmd += 'clname=\"PULSE_NUMBER\" expr=\"%s\" rowrange=\"-\" tform="K"' % operation2
		print(cmd); os.system(cmd)		

		cmd  = 'fcalc clobber=yes infile=%s+1 ' % outfits 
		cmd += 'outfile=%s ' % outfits
		cmd += 'clname=\"PULSE_PHASE\" expr=\"%s\" rowrange=\"-\"' % operation3
		print(cmd); os.system(cmd)	

		cmd  = 'fcalc clobber=yes infile=%s+1 ' % outfits 
		cmd += 'outfile=%s ' % outfits
		cmd += 'clname=\"MOD_PULSE_NUMBER\" expr=\"%s\" rowrange=\"-\" tform="K"' % operation4
		print(cmd); os.system(cmd)				

		#cmd  = 'fcalc clobber=yes infile=%s+1 ' % outfits 
		#cmd += 'outfile=%s ' % outfits
		#cmd += 'clname=\"PULSE_PHASE\" expr=\"%s\" rowrange=\"-\"' % operation2
		#print cmd; os.system(cmd)

		#cmd  = 'fcalc clobber=yes infile=%s+1 ' % outfits 
		#cmd += 'outfile=%s ' % outfits
		#cmd += 'clname=\"PULSE_PHASE\" expr=\"%s\" rowrange=\"-\"' % operation3
		#print cmd; os.system(cmd)

		#cmd  = 'fcalc clobber=yes infile=%s+1 ' % outfits 
		#cmd += 'outfile=%s ' % outfits
		#cmd += 'clname=\"PULSE_NUMBER\" expr=\"%s\" rowrange=\"-\"' % operation4
		#print cmd; os.system(cmd)		

		f = open('temp_header.txt','w')
		f.write(history_dump)
		f.close()
		cmd  = ''
		for i in range(0,2):
			cmd += 'fthedit %s+%d @temp_header.txt\n' % (outfits,i)
			cmd += 'fparkey %.12f "%s[%d]" NU comm="Frequency for the PULSE_PHAES column" add=yes\n' % (nu,outfits,i)
			cmd += 'fparkey %.7e "%s[%d]" NUDOT comm="Frequency derivative for the PULSE_PHAES column" add=yes\n' % (nudot,outfits,i)	
			cmd += 'fparkey %.7e "%s[%d]" NU2DOT comm="Frequency 2nd derivative for the PULSE_PHAES column" add=yes\n' % (nu2dot,outfits,i)	
			cmd += 'fparkey %.7e "%s[%d]" NU3DOT comm="Frequency 3rd derivative for the PULSE_PHAES column" add=yes\n' % (nu3dot,outfits,i)	
			cmd += 'fparkey %.7e "%s[%d]" NU4DOT comm="Frequency 4th derivative for the PULSE_PHAES column" add=yes\n' % (nu4dot,outfits,i)						
			cmd += 'fparkey %.6f "%s[%d]" EPOCH comm="EPOCH for the PULSE_PHAES column" add=yes\n' % (epoch,outfits,i)
		cmd += 'rm -f temp_header.txt'
		print(cmd); os.system(cmd)

if __name__=="__main__":

	usage = """
This python script calcurates pulse phases (PULSE_PHASE) of individual X-ray photons, 
based on an ephemeris of rotation frequency nu (Hz), its derivative (Hz/sec), and the 
folding epoch in a unit of TIME column. For negative optional value, use an equal 
symbol, e.g., --nudot=-3.687105e-10. We employed follogin formula. 
	nu*(t-t0) + nudot*(t-t0)^2/2.0 + nu2dot*(t-t0)^3/6.0 + nu3dot*(t-t0)^4/24.0 + nu4dot*(t-t0)^5/120.0
where t is an event time (TIME column), t0 is an epoch for folding. 

(Note) For folding with period and its derivative
when constant P and Pdot are given, phases when each X-ray photons comes are calcurated
by integrtion of 1/(P + t*Pdot) from 0 to t. Here t is a difference in time of the arrival 
time and the epoch. This gives, 1/(Pdot){ ln(1 + Pdot/P *t) }. 

(Reference); 
* http://heasarc.gsfc.nasa.gov/docs/xte/abc/time.html
* http://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
* http://www.starlink.rl.ac.uk/star/docs/sun67.htx/node226.html
"""

	parser = argparse.ArgumentParser(
		prog='faddphase_nu.py',
		usage='faddphase_nu.py infits epoch nu [--outfits] [--nudot] [--nudot2] [--nudot3] [--nudot3] [--nudot4] [--offset]',
		description=usage,
		epilog='',
		add_help=True)

	parser.add_argument('infits',metavar='infits',type=str,        
		help='Input event fits file.')	
	parser.add_argument('epoch',metavar='epoch',type=float,        
		help='Folding epoch in a unit of TIME.')
	parser.add_argument('nu',metavar='nu',type=float,        
		help='Folding frequency nu (Hz)')
	parser.add_argument('--nudot',metavar='nudot',type=float, default=0.0,  
		help='Folding frequency derivative nudot (Hz/s)')	
	parser.add_argument('--nu2dot',metavar='nu2dot',type=float, default=0.0,  
		help='Folding frequency 2nd derivative nu2dot (Hz/s2)')	
	parser.add_argument('--nu3dot',metavar='nu3dot',type=float, default=0.0,  
		help='Folding frequency 3rd derivative nu3dot (Hz/s3)')	
	parser.add_argument('--nu4dot',metavar='nu4dot',type=float, default=0.0,  
		help='Folding frequency 4th derivative nu4dot (Hz/s4)')	
	parser.add_argument('--outfits',metavar='outfits',type=str,default=None,
		help='Output fits file.')		
	parser.add_argument('--offset',metavar='offset',type=str,default="0.0",
		help='Optional offset to phase.')	
	parser.add_argument('--flag_mjd',action='store_true',dest='flag_mjd',
		default=False,help='Flag MJD (BARY_TIME) calculation.')  													
	args = parser.parse_args()	
	#print(args)

	evtfits = EventFits(args.infits)
	evtfits.faddphase_nu(args.epoch,args.nu,
		nudot=float(args.nudot),nu2dot=float(args.nu2dot),
		nu3dot=float(args.nu3dot),nu4dot=float(args.nu4dot),
		outfits=args.outfits,offset=float(args.offset),
		flag_mjd=args.flag_mjd)
	
