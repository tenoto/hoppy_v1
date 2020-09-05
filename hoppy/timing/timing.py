# -*- coding: utf-8 -*-

import os 
import sys
from datetime import datetime 

import numpy as np

from astropy.io import fits

import matplotlib.pylab as plt 

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

	"""
	def plot_profile(self,outfits=None,phase_colname='PULSE_PHASE',
		binsize=0.05):
		sys.stdout.write("--%s--\n" % sys._getframe().f_code.co_name)

		if outfits == None:
			outfits = '%s.fht' % os.path.splitext(os.path.basename(self.eventfits))[0]

		cmd  = 'fhisto '		
		cmd += '%s ' % self.eventfits
		cmd += '%s ' % outfits
		cmd += '%s ' % phase_colname
		cmd += '%.6f ' % binsize
		cmd += 'lowval=0.0 '
		cmd += 'highval=1.0 '
		cmd += 'outcolx="PULSE_PHASE" '
		cmd += 'outcoly="COUNT" '		
		cmd += 'outcolz="COUNT_ERROR" '				
		print(cmd);os.system(cmd)
	"""

	def get_profile(self,nphase=40,phase_colname='PULSE_PHASE'):
		sys.stdout.write("--%s--\n" % sys._getframe().f_code.co_name)

		hdu = fits.open(self.eventfits) 
		phase_array = hdu['EVENTS'].data[phase_colname]
		nevt = len(phase_array)
		profile_count, profile_binedge, profile_patches = plt.hist(phase_array,nphase,range=[0.0,1.0],histtype='step')
		profile_bincenter = 0.5*(profile_binedge[1:]+profile_binedge[:-1])	
		x = np.array(profile_bincenter)
		y = np.array(profile_count)
		xe = np.full(len(x),0.5/float(nphase))
		ye = np.array(np.sqrt(profile_count))
		ny = y / float(nevt)
		nye = ye / float(nevt)
		return x,xe,y,ye,ny,nye

	def plot_two_profiles(self,x1,x1e,y1,y1e,outpdf,
			label1='',label2='',ylabel='Counts',title='',legend_title='',
			xmin=0.0,xmax=2.0,ymin=None,ymax=None):
		sys.stdout.write("--%s--\n" % sys._getframe().f_code.co_name)

		outdir = os.path.dirname(outpdf)
		if outdir != '' and not os.path.exists(outdir):
			os.makedirs(outdir)

		fontsize = 14

		plt.clf()
		fig, axes = plt.subplots(1,1,figsize=(9.6,5.6))
		plt.errorbar(np.hstack((x1,x1+1.0)),np.tile(y1,2),xerr=np.tile(x1e,2),yerr=np.tile(y1e,2),
			marker='',color='k',drawstyle='steps-mid',linewidth=1.0,label=label1)
		#if ymax != None:
		#	plt.vlines([1.0],0.0,ymax,'k',linestyles='dashed',linewidth=1.0)  
		plt.title(title)
		legend = axes.legend(loc='upper right',shadow=False,fontsize=11.0,title=legend_title)

		axes.set_xlim(float(xmin),float(xmax))
		if ymin != None and ymax != None:
			axes.set_ylim(float(ymin),float(ymax))
		axes.set_xlabel('Pulse Phase',fontsize=fontsize)			
		axes.set_ylabel(ylabel,fontsize=fontsize)
		axes.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
		axes.tick_params(axis="both", which='major', direction='in', length=5)
		axes.tick_params(axis="both", which='minor', direction='in', length=3)				
		plt.subplots_adjust(wspace=0, hspace=0)
		plt.tight_layout(pad=2)		
		plt.rcParams["font.family"] = "serif"
		plt.rcParams["mathtext.fontset"] = "dejavuserif"				
		plt.savefig(outpdf,dpi=300)
	
	def extract_phase_spectrum(self,outpha,filter_expression,exposure_fraction=1.0):
		sys.stdout.write("--%s--\n" % sys._getframe().f_code.co_name)

		cmd  = 'rm -f tmp_fselect.evt xselect.log'
		print(cmd);os.system(cmd)

		cmd  = 'fselect '
		cmd += '%s ' % self.eventfits
		cmd += 'tmp_fselect.evt '
		cmd += '"%s" ' % filter_expression
		print(cmd);os.system(cmd)

		cmd = 'rm -f xselect.log'
		print(cmd);os.system(cmd)

		cmd  = 'xselect<<EOF\n'
		cmd += 'xsel\n'
		cmd += 'read event tmp_fselect.evt ./\n' 
		cmd += 'yes\n'
		cmd += 'extract spec\n'
		cmd += 'save spec %s\n' % outpha
		cmd += 'yes\n'
		cmd += 'exit\n'
		cmd += 'no\n'
		cmd += 'exit\n'	
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)	

		cmd = 'rm -f xselect.log'
		print(cmd);os.system(cmd)

		cmd  = 'rm -f tmp_fselect.evt'
		print(cmd);os.system(cmd)

		hdu = fits.open(outpha)
		exposure = hdu[0].header['EXPOSURE']
		new_exposure = exposure * exposure_fraction
		for i in range(len(hdu)):
			cmd = 'fparkey %.8f %s[%d] EXPOSURE' % (new_exposure,outpha,i)
			print(cmd);os.system(cmd)
		cmd = 'fkeyprint %s EXPOSURE' % outpha
		print(cmd);os.system(cmd)





	










