#!/usr/bin/env python

__name__    = 'fxselect_extract_spec'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 May 8'

import os
import sys
#import pyfits
import astropy.io.fits as pyfits
from optparse import OptionParser
from datetime import datetime 

parser = OptionParser()
parser.add_option("-i","--inputevtfits",dest="inputevtfits",
	action="store",help="input event fits file",type="string")
parser.add_option("-o","--outputpha",dest="outputpha",
	action="store",help="output pha file",type="string")
parser.add_option("-r","--rmf",dest="rmf",default=None,
	action="store",help="rmf file",type="string")
parser.add_option("-a","--arf",dest="arf",default=None,
	action="store",help="arf file",type="string")
parser.add_option("-g","--gtifile",dest="gtifile",default=None,
	action="store",help="gti file",type="string")

parser.add_option("-d","--phasemin",dest="phasemin",default=None,
	action="store",help="phase min",type="float")
parser.add_option("-u","--phasemax",dest="phasemax",default=None,
	action="store",help="phase max",type="float")
(options, args) = parser.parse_args()

if options.inputevtfits == None:
	print "input event fits file is needed. %> -i inputevtfits"
	quit()
if options.outputpha == None:	
	print "output pha file is needed. %> -o outputpha"
	quit()

print "input event fits file  : %s " % options.inputevtfits
print "output pha file: %s " % options.outputpha	

if not os.path.exists(options.inputevtfits):
	print "input event fits file does not exists: %s" % options.inputevtfits
	quit()
if os.path.exists(options.outputpha):
	print "output pha file has already existed: %s " % options.outputpha
	quit()

hdu = pyfits.open(options.inputevtfits)
if len(hdu['EVENTS'].data) == 0:
	sys.stdout.write('Skip: No events in fits file %s' % options.inputevtfits)
	quit()

if options.phasemin != None and options.phasemax != None:
	exposure = float(hdu['EVENTS'].header['EXPOSURE'])
	new_exposure = (options.phasemax - options.phasemin)*exposure
	cmd  = 'rm -f tmp_sel.evt\n'
	cmd += 'fselect %s tmp_sel.evt expr="(PULSE_PHASE>=%.5f)&&(PULSE_PHASE<=%.5f)"' % (options.inputevtfits,options.phasemin,options.phasemax)
	print(cmd);os.system(cmd)
	cmd = """
xselect <<EOF
xsel
read event tmp_sel.evt .
yes
set phaname PI
"""
	if options.gtifile != None:
		cmd += "filter time file %s\n" % options.gtifile
		cmd += "show filter\n"
	cmd += """extract spectrum
save spectrum
%s
exit
no
EOF
""" % (options.outputpha)
	print(cmd);os.system(cmd)	
	cmd  = "fparkey %.5f %s[0] EXPOSURE\n" % (new_exposure,options.outputpha)
	cmd += "fparkey %.5f %s[1] EXPOSURE\n" % (new_exposure,options.outputpha)	
	cmd += "fparkey %.5f %s[2] EXPOSURE\n" % (new_exposure,options.outputpha)		
	print(cmd);os.system(cmd)
else:
	cmd = """
xselect <<EOF
xsel
read event %s .
yes
set phaname PI
""" % options.inputevtfits
	if options.gtifile != None:
		cmd += "filter time file %s\n" % options.gtifile
		cmd += "show filter\n"
	cmd += """extract spectrum
save spectrum
%s
exit
no
EOF
""" % (options.outputpha)
	print(cmd);os.system(cmd)

cmd = 'rm -f xselect.log xsel_timefile.asc'
print(cmd);os.system(cmd)

if options.rmf != None:
	cmd = 'fparkey %s %s RESPFILE add=yes' % (options.rmf, options.outputpha)
	print(cmd);os.system(cmd)
if options.arf!= None:
	cmd = 'fparkey %s %s ANCRFILE add=yes' % (options.arf, options.outputpha)	
	print(cmd);os.system(cmd)
