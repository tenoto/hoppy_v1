#!/usr/bin/env python

__name__    = 'fplot_pulseprofile'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 March 30'

import os
import sys
#import pyfits
import astropy.io.fits as pyfits
import subprocess
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i","--inputevtfits",dest="inputevtfits",
	action="store",help="input event fits file, PHASE column is required",type="string")
parser.add_option("-o","--outputhistfits",dest="outputhistfits",
	action="store",help="output histogram fits file",type="string")
parser.add_option("-n","--nbin",dest="nbin",
	action="store",help="Phase bin for a one pulse",type="int")
parser.add_option("-c","--column",dest="column",default="PULSE_PHASE",
	action="store",help="Column name for PHASE (Defaul=PULSE_PHASE)",type="string")
parser.add_option("-t","--title",dest="title",default="",
	action="store",help="title",type="string")
parser.add_option("-b","--backgroundrate",dest="backgroundrate",default=0.0,
	action="store",help="Background rate (cps)",type="float")
(options, args) = parser.parse_args()

if options.inputevtfits == None:
	print "input event fits file is needed. %> -i inputevtfits"
	quit()
if options.outputhistfits == None:	
	print "output histogram fits file is needed. %> -o outputhistfits"
	quit()
if options.nbin == None:	
	print "Phase bin for a one pulse is needed. %> -n nbin"
	quit()

print "input event fits file  : %s " % options.inputevtfits
print "output event fits file : %s " % options.outputhistfits	
print "nbin : %.2f " % options.nbin

if not os.path.exists(options.inputevtfits):
	print "input event fits file does not exists: %s" % options.inputevtfits
	quit()
if os.path.exists(options.outputhistfits):
	print "output histogram fits file has already existed: %s " % options.outputhistfits
	quit()

basename = os.path.splitext(options.outputhistfits)[0]

tmp1 = '%s_tmp1.fht' % basename
cmd  = 'rm -f %s' % tmp1
print(cmd);os.system(cmd)

phase_bin = 1.0/float(options.nbin)
cmd  = 'fhisto %s %s ' % (options.inputevtfits,tmp1)
cmd += '%s %.3f ' % (options.column, phase_bin)
cmd += 'lowval=0.0 highval=1.0 outcolx=%s outcoly=COUNTS outcolz=ERROR extname=PROFILE' % options.column
print(cmd);os.system(cmd)

tmp2 = '%s_tmp2.fht' % basename
cmd  = 'fcalc clobber=yes %s %s %s %s+1.0' % (tmp1,tmp2,options.column,options.column)
print(cmd);os.system(cmd)

tmp3 = '%s_tmp3.fht' % basename
cmd = 'ftmerge %s,%s %s' % (tmp1,tmp2,tmp3)
print(cmd);os.system(cmd)

cmd = 'get_total_gti_exposure.py %s' % options.inputevtfits
dump = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
exposure = float(dump.stdout.readlines()[0].strip())
cmd  = 'fparkey %.5f %s[1] EXPOSURE add=yes\n' % (exposure,tmp3)
cmd += 'fcalc %s %s RATE "COUNTS/(#EXPOSURE/((float) %d ))"\n' % (tmp3,options.outputhistfits,options.nbin)
#cmd += 'fcalc clobber=yes %s %s RATE_ERR "sqrt(COUNTS)/(#EXPOSURE/((float) %d ))" ' % (tmp3,options.outputhistfits,options.nbin)
cmd += 'fcalc %s %s RATE_ERR "sqrt(COUNTS)/(#EXPOSURE/((float) %d ))" clobber=YES;' % (options.outputhistfits,options.outputhistfits,options.nbin)
print(cmd);os.system(cmd)

hdu = pyfits.open(options.outputhistfits)
rate_list = hdu['PROFILE'].data['RATE'][0:options.nbin]
max_rate = max(rate_list)
min_rate = min(rate_list)
bgdsub_max_rate = max_rate - options.backgroundrate
bgdsub_min_rate = min_rate - options.backgroundrate
pulsed_fraction = (bgdsub_max_rate-bgdsub_min_rate)/(bgdsub_max_rate+bgdsub_min_rate)
cmd  = 'fparkey %.5f %s[1] MAX_WBGD add=yes comm="Maximum rate (cps) before background subtraction.";\n' % (max_rate,options.outputhistfits)
cmd += 'fparkey %.5f %s[1] MIN_WBGD add=yes comm="Minimum rate (cps) before background subtraction.";\n' % (min_rate,options.outputhistfits)
cmd += 'fparkey %.5f %s[1] MAX_RATE add=yes comm="Background-subtracted maximum rate (cps).";\n' % (bgdsub_max_rate,options.outputhistfits)
cmd += 'fparkey %.5f %s[1] MIN_RATE add=yes comm="Background-subtracted minimum rate (cps).";\n' % (bgdsub_min_rate,options.outputhistfits)
cmd += 'fparkey %.5f %s[1] PULFRACT add=yes comm="Pulsed fraction after the background subtraction.";\n' % (pulsed_fraction,options.outputhistfits)
cmd += 'fparkey %d %s[1] NBINPHAS add=yes comm="Number of phase bins";\n' % (options.nbin,options.outputhistfits)
cmd += 'fparkey %.5e %s[1] BACKRATE add=yes comm="Subtracted background rate";\n' % (options.backgroundrate,options.outputhistfits)
print(cmd);os.system(cmd)
subtitle = 'BGD=%.2e (cps) Max=%.2e Min=%.2e Pulsed_Fraction=%.1f%%' % (options.backgroundrate,
	bgdsub_max_rate,bgdsub_min_rate,pulsed_fraction*100.0)

psfile = '%s.ps' % os.path.splitext(options.outputhistfits)[0]
cmd  = 'fplot %s %s "RATE[RATE_ERR]" - /xw @ ' % (options.outputhistfits,options.column)
cmd += '<<EOF\n'
cmd += 'line step on\n'
cmd += 'time off\n'
cmd += 'lwid 5 on\n'
cmd += 'lwid 5\n'
cmd += 'la t %s\n' % options.title
cmd += 'la f %s\n' % subtitle 
cmd += 'hard %s/cps\n' % psfile
cmd += 'quit\n'
cmd += 'EOF\n'
print(cmd);os.system(cmd)

pdffile = '%s.pdf' % os.path.splitext(os.path.basename(psfile))[0]
outdir = os.path.dirname(options.outputhistfits)
cmd  = 'ps2pdf %s ;' % (psfile)
if outdir == '':
	cmd += 'mv %s . ' % pdffile
else:	
	cmd += 'mv %s %s ' % (pdffile,outdir)
print(cmd);os.system(cmd)

cmd = 'rm -f %s %s %s %s' % (tmp1,tmp2,tmp3,psfile)
print(cmd);os.system(cmd)



