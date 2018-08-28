#!/usr/bin/env python

import sys 
import math 

# see p.37
# https://www.dropbox.com/s/xlacg6vaz31ft9w/AHCookBook_v1.02.pdf?dl=0

if len(sys.argv) != 3 :
	sys.stderr.write('%s enegry_keV velocity_kms\n' % sys.argv[0])
	exit()

velocity_kms_light = 299792.458 # km/s
energy_keV = float(sys.argv[1])
velocity_kms = float(sys.argv[2])
theta = 0.0

beta = velocity_kms / velocity_kms_light
fraction = (1 - beta * math.cos(theta))/math.sqrt(1-beta**2) 
shifted_energy = fraction * energy_keV
print("fraction:%.6f" % fraction)
print("shifted energy: %.6f keV" % shifted_energy)
