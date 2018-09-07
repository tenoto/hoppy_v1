#!/usr/bin/env python

import sys 
import astropy.io.fits as pyfits

if len(sys.argv) != 2:
	sys.stderr.write("usage: %s infits\n" % sys.argv[0])
	quit()
infits = sys.argv[1]

hdu = pyfits.open(infits)
print(len(hdu['EVENTS'].data))
