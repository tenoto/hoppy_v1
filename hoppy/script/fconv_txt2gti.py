#!/usr/bin/env python

__name__    = 'fconv_txt2gti'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 June 7'

import os
import sys 
import astropy.io.fits as pyfits 
from optparse import OptionParser
from datetime import datetime 

parser = OptionParser()
parser.add_option("-i","--inputfile",dest="inputfile",
	action="store",help="input text file, tab separated",type="string")
parser.add_option("-o","--outputgti",dest="outputgti",
	action="store",help="output gti file name",type="string")
(options, args) = parser.parse_args()

if options.inputfile == None:
	sys.stderr.write("input text file is needed. %> -i inputfile")
	quit()
if options.outputgti == None:	
	sys.stderr.write("output gti file name is needed. %> -o outputgti")
	quit()
sys.stdout.write("inputfile : %s " % options.inputfile)
sys.stdout.write("outputgti: %s " % options.outputgti)

if not os.path.exists(options.inputfile):
	sys.stderr.write("input text fits file does not exists: %s" % options.inputfile)
	quit()
if os.path.exists(options.outputgti):
	sys.stderr.write("output gti file has already existed: %s " % options.outputgti)
	quit()


f = open('gti_columns.txt','w')
dump = """START D s
STOP D s

"""
f.write(dump)
f.close()

f = open('gti_header.txt','w')
dump = """XTENSION= 'BINTABLE'           / binary table extension
BITPIX  =                    8 / 8-bit bytes
NAXIS   =                    2 / 2-dimensional binary table
NAXIS1  =                   16 / width of table in bytes
NAXIS2  =                    1 / number of rows in table
PCOUNT  =                    0 / size of special data area
GCOUNT  =                    1 / one data group (required keyword)
TFIELDS =                    2 / number of fields in each row
TTYPE1  = 'START   '           / lower GTI boundary
TFORM1  = 'D       '           / data format of field: 8-byte DOUBLE
TUNIT1  = 's       '           / physical unit of field
TTYPE2  = 'STOP    '           / upper GTI boundary
TFORM2  = 'D       '           / data format of field: 8-byte DOUBLE
TUNIT2  = 's       '           / physical unit of field
EXTNAME = 'STDGTI  '           / The name of this table
HDUCLASS= 'OGIP    '           / format conforms to OGIP standard
HDUCLAS1= 'GTI     '           / table contains Good Time Intervals
HDUCLAS2= 'STANDARD'           / standard Good Time Interval table
ONTIME  = 0.00000000000000E+00 / [s] sum of all Good Time Intervals
TSTART  = 0.00000000000000E+00 / [s] Lower bound of first GTI
TSTOP   = 0.00000000000000E+00 / [s] Uppler bound of last GTI
TIMEUNIT= 's       '           / All times in s unless specified otherwise
TIMESYS = 'TT      '           / XMM time will be TT (Terrestial Time)
TIMEREF = 'LOCAL   '           / Reference location of photon arrival times
TASSIGN = 'SATELLITE'          / Location of time assignment
TIMEZERO=                    0 / Clock correction (if not zero)
CLOCKAPP=                    T / Clock correction applied?
MJDREFI =                56658 / MJD reference day
MJDREFF = 7.775925925925930E-04 / MJD reference (fraction of day)
"""
f.write(dump)
f.close()


f = open('tmp_gti_data_shrink.txt','w')
flag_first = True
for line in open(options.inputfile):
	cols = line.split()
	if cols[0] == '#':
		continue
	tmp_TSTART = cols[0]
	tmp_TSTOP  = cols[1]
	if flag_first:
		TSTART = tmp_TSTART
		prev_TSTOP  = tmp_TSTOP
		flag_first = False
		continue

	if tmp_TSTART != prev_TSTOP:
		dump = '%s %s\n' % (TSTART,prev_TSTOP)
		f.write(dump)
		TSTART = tmp_TSTART 
	prev_TSTOP = tmp_TSTOP
dump = '%s %s\n' % (TSTART,prev_TSTOP)
f.write(dump)
f.close()

#cmd  = 'ftcreate gti_columns.txt %s %s headfile=gti_header.txt extname="GTI" clobber=yes\n' % (options.inputfile,options.outputgti)
cmd  = 'ftcreate gti_columns.txt tmp_gti_data_shrink.txt %s headfile=gti_header.txt extname="GTI" clobber=yes\n' % (options.outputgti)
print(cmd);os.system(cmd)

cmd  = 'rm -f tmp_gti_data_shrink.txt gti_columns.txt gti_header.txt'
print(cmd);os.system(cmd)


