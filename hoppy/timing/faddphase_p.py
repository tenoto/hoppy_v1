#!/usr/bin/env python
# 
# This script add "PHASE" column calculated from "PERIOD" column.
# The original script is aeplsphase.pl created by T.Enoto 
# created by Teruaki Enoto 2015-02-07 
# converted from perl script 

__name__    = 'faddphase'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2015 Feb. 7'

import os
from optparse import OptionParser
from datetime import datetime 

parser = OptionParser()
parser.add_option("-i","--inputfits",dest="inputfits",
	action="store",help="input fits file",type="string")
parser.add_option("-o","--outputfits",dest="outputfits",
	action="store",help="output fits file",type="string")
parser.add_option("-p","--period",dest="period",
	action="store",help="period (s)",type="float")
parser.add_option("-e","--epoch",dest="epoch",
	action="store",help="epoch (s)",type="float")
parser.add_option("-d","--pdot",dest="pdot",
	action="store",help="period derivative (s/s)",type="float", default=0.0)
(options, args) = parser.parse_args()

if options.inputfits == None:
	print "input fits file is needed. %> -i inputfits"
	quit()
if options.outputfits == None:	
	print "output fits file is needed. %> -o outputfits"
	quit()
if options.period == None:
	print "period (s) is needed. %> -p period"
	quit()
if options.epoch == None:
	print "epoch (s) is needed. %> -e epoch"
	quit()
print "inputfits : %s " % options.inputfits
print "outputfits: %s " % options.outputfits
print "period: %.12f (s)" % options.period
print "pdot: %.8e (s/s)" % options.pdot
print "epoch: %.6f (s)" % options.epoch

if not os.path.exists(options.inputfits):
	print "input file does not exists: %s" % options.inputfits
	quit()
if os.path.exists(options.outputfits):
	print "output file has already existed: %s " % options.outputfits
	quit()

if options.pdot == None:
	options.pdot = 0.0

if options.pdot == 0:
	operation = "((TIME-%.6f) %% %.12f) / %.12f" % (
		options.epoch, options.period, options.period)
else:
	operation = "(log(1.0+%.8e*(TIME-%.6f)/%.12f))/%.8e %% 1.0" % (
		options.pdot, options.epoch, options.period, options.pdot)

# ----- Reference of the calcuration -----
# When constant P and Pdot are given, phases when each X-ray photons
# comes are calcurated by integrtion of 1/(P + t*Pdot) from 0 to t.
# Here t is a difference in time of the arrival time and the epoch.
# This gives,
#     1/(Pdot){ ln(1 + Pdot/P *t) }
# ----------------------------------------
print operation
cmd  = 'fcalc infile=%s+1 ' % options.inputfits
cmd += 'outfile=%s ' % options.outputfits
cmd += 'clname=\"PHASE\" expr=\"%s\" rowrange=\"-\"' % operation
print cmd; os.system(cmd)

operation2 = ' PHASE<0 ? PHASE+1 : PHASE '
cmd  = 'fcalc clobber=yes infile=%s ' % options.outputfits
cmd += 'outfile=%s ' % options.outputfits
cmd += 'clname=\"PHASE\" expr=\"%s\" rowrange=\"-\"' % operation2
print cmd; os.system(cmd)

out = """
HISTORY -----------------------------------------------------
HISTORY  %s version %s at %s
HISTORY -----------------------------------------------------
HISTORY   inputfits='%s'
HISTORY   outputfits='%s'
HISTORY   period = %.12f / Period (sec) used at phase calcuration
HISTORY   pdot  = %.8e   / Pdot (sec/sec) used at phase calcuration
HISTORY   epoch = %.6f   / Epoch (s) used at phase calcuration
HISTORY   %s 
HISTORY   %s 
""" % (__name__, __version__, datetime.now().strftime('%Y-%m-%dT%H:%M:%S'),
	options.inputfits, options.outputfits,
	options.period, options.pdot, options.epoch,
	operation, operation2)
print out
f = open('temp_header.txt','w')
f.write(out)
f.close()
cmd  = ''
for i in range(0,2):
	cmd += 'fthedit %s+%d \@temp_header.txt\n' % (options.outputfits,i)
	cmd += 'fparkey %.12f "%s[%d]" PERIOD comm="Period for the PHASE column" add=yes\n' % (
		options.period, options.outputfits, i)
	cmd += 'fparkey %.8e "%s[%d]" PDOT comm="Period derivative for the PHASE column" add=yes\n' % (
		options.pdot, options.outputfits, i)	
	cmd += 'fparkey %.6f "%s[%d]" EPOCH comm="EPOCH for the PHASE column" add=yes\n' % (
		options.epoch, options.outputfits, i)	
cmd += 'rm -f temp_header.txt'
print cmd; os.system(cmd)

"""
NAME
        aeplsphase (${ver})

DESCRIPTION
        This perl script calcurates pulse phases of individual X-ray
        photons, based on the period (sec), its derivative (ses/sec)
        and the epoch (MJD in TT; Terrestrial Time). This script calls
        aetimecalc and convert MJD to Suzaku time. Used parameters 
	are written to the header file. 

        When constant P and Pdot are given, phases when each X-ray photons
	comes are calcurated by integrtion of 1/(P + t*Pdot) from 0 to t.
	Here t is a difference in time of the arrival time and the epoch.
	This gives,
	     1/(Pdot){ ln(1 + Pdot/P *t) }
        We employed the following conversion fomula in the same way as XRONOS.
	  Pdot  = 0 ==> "((TIME-\$epoch) \% \$per) / \$per"
	  Pdot != 0 ==> "(log(1.0+\$pdot*(TIME-\$epoch)/\$per))/\$pdot \% 1.0"
	Reference; 
            http://heasarc.gsfc.nasa.gov/docs/xte/abc/time.html
            http://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
            http://www.starlink.rl.ac.uk/star/docs/sun67.htx/node226.html

AUTHOR
        Teruaki Enoto (enoto\@juno.phys.s.u-tokyo.ac.jp)

SYNTAX 
        $0 -i infits -o outfits -p period(s) -dp dpdot(ss-1) -e epoch

OPTIONS
        -p      Period (s)

	-dp     Period derivative (ss\-1)

	-e      Epoch (MJD in TT; Modified Julian Datein of Terrestrial Time)

        -v      Verbosity: print detailed informations during the 
                process.
               
        -h      Give syntax. This argument must be the first.

        -help   Give detailed help (print this message). 
		This argument must be the first.

EXAMPLE
        .\/aeplsphase \\
	   -i crab.evt -o crab_phase.evt \\
	   -p 3.358087649664188e-02 -dp 4.2059137835e-13 \\
	   -e 53628.000000367673611111
EndHelp
"""
