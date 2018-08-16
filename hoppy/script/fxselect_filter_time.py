#!/usr/bin/env python

__name__    = 'fxselect_filter_time'
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
parser.add_option("-o","--outputevtfits",dest="outputevtfits",
	action="store",help="output lc file",type="string")
parser.add_option("-d","--tstart",dest="tstart",
	action="store",help="start time (Spacecraft time)",type="float")
parser.add_option("-u","--tstop",dest="tstop",
	action="store",help="stop time (Spacecraft time)",type="float")
(options, args) = parser.parse_args()

if options.inputevtfits == None:
	print "input event fits file is needed. %> -i inputevtfits"
	quit()
if options.outputevtfits == None:	
	print "output lc file is needed. %> -o outputevtfits"
	quit()
if options.tstart == None:	
	print "tstart is needed. %> -o tstart"
	quit()
if options.tstop == None:	
	print "tstop is needed. %> -o tstop"
	quit()

print "input event fits file  : %s " % options.inputevtfits
print "output event fits file : %s " % options.outputevtfits	

if not os.path.exists(options.inputevtfits):
	print "input event fits file does not exists: %s" % options.inputevtfits
	quit()
if os.path.exists(options.outputevtfits):
	print "output event fits file has already existed: %s " % options.outputevtfits
	quit()

hdu = pyfits.open(options.inputevtfits)
if len(hdu['EVENTS'].data) == 0:
	sys.stdout.write('Skip: No events in fits file %s' % options.inputevtfits)
	quit()

cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc
%.8f,%.8f
x
show filter 
extract event
save event 
%s
no
exit
no
EOF
""" % (options.inputevtfits,options.tstart,options.tstop,options.outputevtfits)
print(cmd);os.system(cmd)	

cmd = 'rm -f xselect.log xsel_timefile.asc'
print(cmd);os.system(cmd)

