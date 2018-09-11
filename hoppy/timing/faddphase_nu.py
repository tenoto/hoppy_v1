#!/usr/bin/env python
# 
# This script add "PHASE" column calculated from nu, nudot, nuddot, and epoch. 
# The original script is aeplsphase.pl created by T.Enoto 
# created by Teruaki Enoto 2015-02-07 
# converted from perl script 

__name__    = 'faddphase_nu'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2017 Aug. 26'

import os
from optparse import OptionParser
from datetime import datetime 

parser = OptionParser()
parser.add_option("-i","--inputfits",dest="inputfits",
	action="store",help="input fits file",type="string")
parser.add_option("-o","--outputfits",dest="outputfits",
	action="store",help="output fits file",type="string")
parser.add_option("-e","--epoch",dest="epoch",
	action="store",help="epoch (s)",type="float")
parser.add_option("-n","--nu",dest="nu",
	action="store",help="nu (Hz)",type="float")
parser.add_option("-s","--nudot",dest="nudot",
	action="store",help="nu derivative (Hz/s)",type="float", default=0.0)
parser.add_option("-d","--nudot2",dest="nudot2",
	action="store",help="nu 2nd derivative (Hz/s2)",type="float", default=0.0)
parser.add_option("-t","--nudot3",dest="nudot3",
	action="store",help="nu 3rd derivative (Hz/s3)",type="float", default=0.0)
parser.add_option("-q","--nudot4",dest="nudot4",
	action="store",help="nu 4th derivative (Hz/s4)",type="float", default=0.0)
parser.add_option("-a","--offset",dest="offset",
	action="store",help="offset (phase)",type="float", default=0.0)
(options, args) = parser.parse_args()

if options.inputfits == None:
	print "input fits file is needed. %> -i inputfits"
	quit()
if options.outputfits == None:	
	print "output fits file is needed. %> -o outputfits"
	quit()
if options.nu == None:
	print "nu (Hz) is needed. %> -n nu"
	quit()
if options.epoch == None:
	print "epoch (s) is needed. %> -e epoch"
	quit()
print "inputfits : %s " % options.inputfits
print "outputfits: %s " % options.outputfits
print "nu: %.12e (s)" % options.nu
print "epoch: %.7f (s)" % options.epoch
print "nudot: %.7e (s)" % options.nudot
print "nudot2: %.7e (s)" % options.nudot2
print "nudot3: %.7e (s)" % options.nudot3
print "nudot4: %.7e (s)" % options.nudot4
print "offset: %.7f (s)" % options.offset

if not os.path.exists(options.inputfits):
	print "input file does not exists: %s" % options.inputfits
	quit()
if os.path.exists(options.outputfits):
	print "output file has already existed: %s " % options.outputfits
	quit()	

operation = "%.7f + (%.12e*(TIME-(%.7f)) + (%.7e)*(TIME-(%.7f))**2/2.0 + (%.7e)*(TIME-(%.7f))**3/6.0 + (%.7e)*(TIME-(%.7f))**4/24.0 + (%.7e)*(TIME-(%.7f))**5/120.0) %% 1.0" % (
	options.offset, options.nu, options.epoch, options.nudot, options.epoch, options.nudot2, options.epoch, options.nudot3, options.epoch, options.nudot4, options.epoch) 
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
HISTORY   nu = %.12f / Frequency (Hz) used at phase calcuration
HISTORY   nudot  = %.7e   / Frequency derivative (Hz/sec) used at phase calcuration
HISTORY   nudot2 = %.7e   / Frequency 2nd derivative (Hz/sec2) used at phase calcuration
HISTORY   nudot3 = %.7e   / Frequency 3rd derivative (Hz/sec3) used at phase calcuration
HISTORY   nudot4 = %.7e   / Frequency 4th derivative (Hz/sec4) used at phase calcuration
HISTORY   epoch  = %.7f   / poch (s) used at phase calcuration
HISTORY   offset = %.7e   / Phase offset 
HISTORY   %s 
HISTORY   %s 
""" % (__name__, __version__, datetime.now().strftime('%Y-%m-%dT%H:%M:%S'),
	options.inputfits, options.outputfits,
	options.nu, options.nudot, options.nudot2, options.nudot3, options.nudot4, 
	options.epoch, options.offset,
	operation, operation2)
print out
f = open('temp_header.txt','w')
f.write(out)
f.close()
cmd  = ''
for i in range(0,2):
	cmd += 'fthedit %s+%d \@temp_header.txt\n' % (options.outputfits,i)
	cmd += 'fparkey %.12f "%s[%d]" NU comm="Frequency for the PHASE column" add=yes\n' % (
		options.nu, options.outputfits, i)
	cmd += 'fparkey %.7e "%s[%d]" NUDOT comm="Frequency derivative for the PHASE column" add=yes\n' % (
		options.nudot, options.outputfits, i)	
	cmd += 'fparkey %.7e "%s[%d]" NUDOT2 comm="Frequency 2nd derivative for the PHASE column" add=yes\n' % (
		options.nudot2, options.outputfits, i)	
	cmd += 'fparkey %.7e "%s[%d]" NUDOT3 comm="Frequency 3rd derivative for the PHASE column" add=yes\n' % (
		options.nudot3, options.outputfits, i)	
	cmd += 'fparkey %.7e "%s[%d]" NUDOT4 comm="Frequency 4th derivative for the PHASE column" add=yes\n' % (
		options.nudot4, options.outputfits, i)						
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