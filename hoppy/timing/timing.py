# -*- coding: utf-8 -*-

import os 
import sys
from datetime import datetime 

__version__ = '2.01'

class EventFits():
	def __init__(self,eventfits):
		sys.stdout.write("An EventFits object of the HOPPY timing module is generated.\n")
		self.eventfits = eventfits
		sys.stdout.write('Input event fits file: %s\n' % self.eventfits)
		if not os.path.exists(self.eventfits):
			sys.stderr.write('Error: input file does not exists: %s\n' % self.eventfits)

	def faddphase_nu(self,epoch,nu,
		nudot=0.0,nu2dot=0.0,nu3dot=0.0,nu4dot=0.0,
		outfits=None,offset=0.0,flag_mjd=False):
		sys.stdout.write("--%s--\n" % sys._getframe().f_code.co_name)

		if outfits != None:
			outfits = outfits 			
		else:
			outfits = self.eventfits.replace('.evt','_phase.evt').replace('.fits','_phase.fits')
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
		self.eventfits, outfits,
		nu, nudot, nu2dot, nu3dot, nu4dot, epoch, offset, 
		operation1,operation2,operation3)
		print(history_dump)

		cmd  = 'fcalc infile=%s+1 ' % self.eventfits
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
