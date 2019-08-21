#!/usr/bin/env python

import sys
import numpy as np
import astropy.io.fits as pyfits

if len(sys.argv) != 4:
	print("error %s file.arf emin emax" % sys.argv[0])
	quit()

arffile = sys.argv[1]	
emin = float(sys.argv[2])
emax = float(sys.argv[3])

hdu = pyfits.open(arffile)
data = hdu["SPECRESP"].data

flag_emin = (emin <= data["ENERG_LO"])
flag_emax = (emax > data["ENERG_HI"])
flag = np.logical_and(flag_emin,flag_emax)

sum_specresp = sum(data[flag]["SPECRESP"])
num_of_point = sum(flag)

ave_effective_area = float(sum_specresp)/float(num_of_point)
print(ave_effective_area)
