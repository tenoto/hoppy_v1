#!/usr/bin/env python

__name__    = 'xselect_gtifilter'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2020 July 31'

import os
import sys
from optparse import OptionParser

description = """This code is a rapper to extract an event fitsfile with a gti filtering via xselect with an input fits-format gti file.
"""

parser = OptionParser(description=description)
parser.add_option("-i","--inputevtfile",dest="inputevtfile",
	action="store",help="input event file",type="string")
parser.add_option("-g","--inputgtifile",dest="inputgtifile",
	action="store",help="input gti file",type="string")
parser.add_option("-o","--outputevtfile",dest="outputevtfile",
	action="store",help="output event file",type="string")
(options, args) = parser.parse_args()

if options.inputevtfile == None:
	sys.stderr.write("input event file is needed. %> -i inputevtfile")
	quit()
if options.inputgtifile == None:	
	sys.stderr.write("input gti file is needed. %> -g inputgtifile")
	quit()
if options.outputevtfile == None:	
	sys.stderr.write("output event file is needed. %> -g outputevtfile")
	quit()	
sys.stdout.write("inputevtfile : %s " % options.inputevtfile)
sys.stdout.write("inputgtifile: %s " % options.inputgtifile)
sys.stdout.write("outputevtfile: %s " % options.outputevtfile)

if not os.path.exists(options.inputevtfile):
	sys.stderr.write("input event fits file does not exists: %s" % options.inputevtfile)
	quit()
if not os.path.exists(options.inputgtifile):
	sys.stderr.write("input gti fits file does not exists: %s" % options.inputgtifile)
	quit()	
if os.path.exists(options.outputevtfile):
	sys.stderr.write("output event file has already existed: %s " % options.outputevtfile)
	quit()


cmd = 'rm -f xselect.log'
print(cmd);os.system(cmd)

cmd  = 'xselect<<EOF\n'
cmd += 'xsel\n'
cmd += 'read event %s ./\n' % options.inputevtfile
cmd += 'yes\n'
cmd += 'filter time file\n'
cmd += '%s\n' % options.inputgtifile
cmd += 'extract event\n'
cmd += 'save event %s\n' % options.outputevtfile
cmd += 'yes\n'
cmd += 'exit\n'
cmd += 'no\n'
cmd += 'exit\n'	
cmd += 'EOF\n'
print(cmd);os.system(cmd)	

cmd = 'rm -f tmp_gti.txt xselect.log xsel_timefile.asc'
print(cmd);os.system(cmd)

