#!/usr/bin/env python
# Create grppha spectral file (FTOOLS "grppha" command) 
# consistent with the XSPEC rebin command.
# 
# ------------------------------------------
# grppha from XS{EC} "setplot rebin" command
# usage :
#    xspec> setplot rebin 3 10
#    xspec> setplot channel
#    xspec> iplot data
#    PLT> we binned_filename.qdp
#    %> fgrppha.py original.pha binned_filename.qdp [systematic_error]
# output filename : binned_filename.pha
# systematic error (optional): 0.01 means 1% systematic error
# for inclusinon of the systematic error,
# see http://cosmic.riken.jp/suzaku/help/guide/fstep_web/node10.html
# (reproduced from the perl script, 2015-02-27, TE)

import sys
import os 
import math

__author__  = "Teru Enoto"
__version__ = "1.0.0"
__email__   = "teruaki.enoto@nasa.gov"
__status__  = "Fixed"

if len(sys.argv) != 3 and len(sys.argv) != 4:
	sys.stderr.write("usage %> ./fgrppha.py original.pha binned_filename.qdp [systematic_error]\n")
	quit()
if len(sys.argv) == 4:
	flagSystematic = True
	try:
		syserr = float(sys.argv[3])
		sys.stderr.write("Systematic error: %.3f%%\n" % syserr)
	except:
		sys.stderr.write("systematic error should be value: %s\n" % sys.argv[3])
		quit()
else:
	flagSystematic = False	

inpha = sys.argv[1]
inqdp = sys.argv[2]
if not os.path.exists(inpha):
	sys.stderr.write("input pha file %s does not exists.\n" % inpha)
	quit()
if not os.path.exists(inqdp):
	sys.stderr.write("input qdp file %s does not exists.\n" % inqdp)
	quit()
if os.path.splitext(inqdp)[-1] not in ['.qdp']:
	sys.stderr.write("input qdp file must be qdp file (with .qdp extension): %s\n" % inqdp) 
	quit()
outpha = '%s.pha' % os.path.splitext(inqdp)[0]
if os.path.exists(outpha):
	sys.stderr.write("output pha file has already existed: %s\n" % outpha)
	quit()

cmd  = 'grppha<<EOF\n'
cmd += '%s\n' % inpha
cmd += '%s\n' % outpha

flagBody = False
for line in open(inqdp):
	cols = line.split()
	if cols[0] == '!':
		flagBody = True
		continue
	if not flagBody:
		continue
	start   = int(float(cols[0]) - float(cols[1]))
	stop    = int(float(cols[0]) + float(cols[1]) - 1)
	binsize = int(2.0*float(cols[1]))
	cmd += 'group %d %d %d\n' % (start, stop, binsize)
	if flagSystematic:
		binsyserr = math.sqrt(float(binsize)) * syserr
		cmd += 'systematics %d-%d %.4f \n' % (start,stop,binsyserr)

cmd += "exit\n"	
cmd += "EOF\n"	
os.system(cmd)

f = open('temp_header.txt','w')
f.write('HISTORY ---------------\n')
if flagSystematic:
	f.write('HISTORY fgrppha.py %s %s %s\n' % (sys.argv[1],sys.argv[2],sys.argv[3]))
else:	
	f.write('HISTORY fgrppha.py %s %s\n' % (sys.argv[1],sys.argv[2]))
for line in cmd.split('\n'):
	dump = 'HISTORY %s\n' % line
	f.write(dump)
f.close()

cmd = ''
for i in range(2):
	cmd += 'fthedit %s+%d \@temp_header.txt\n' % (outpha,i)
	if flagSystematic:
		cmd += 'fparkey %.4f "%s[%d]" SYSERR comm="systematic error for the spectral bins." add=yes\n' % (
			syserr, outpha, i)
cmd += 'rm -f temp_header.txt'	
os.system(cmd)	
