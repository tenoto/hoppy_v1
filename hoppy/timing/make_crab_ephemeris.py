#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 Sep. 7'
__version__ = '1.00'

code_description = """
This script writes Tempo2-format PAR file from 
the Jodrell Bank Observatory (JBO) Crab monthly ephemeris parameters. 
Original code from Paul Ray. 
"""

import os 
import numpy as np
import argparse 
import astropy.units as u
from astropy.table import Table
from astropy import log

JBO_CRAB_PULSAR_MONTHLY_EPHEMERIS = "http://www.jb.man.ac.uk/pulsar/crab/crab2.txt"

class Ephemeris_CrabPulsar():
	def __init__(self,reference="JBO"):
		if reference == "JBO":
			self.read_JBO_Crab_ephemeris()

	def read_JBO_Crab_ephemeris(self):
		cmd = 'rm -f crab2.txt'
		print(cmd);os.system(cmd)

		cmd = "curl -O %s" % JBO_CRAB_PULSAR_MONTHLY_EPHEMERIS
		print(cmd);os.system(cmd)

		validmonths = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'JLY', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

		day = []; mon = []; yr = []; mjd = []; 
		tJPL = []; tacc = []; 
		f0 = []; f1 = []
		dm = []; dm1 = []
		for line in open('crab2.txt'):
		    line = line.strip()
		    if len(line)==0:
		        continue
		    line = (line.replace('(',' ')).replace(')',' ')
		    cols = line.split()
		    # Make sure second column is a valid month, or skip line. This skips headers, etc.
		    if not cols[1] in validmonths:
		        continue
		    # Skip entries before 55910 when DM1 column was added.
		    if float(cols[3]) < 55910:
		        continue
		    #print(line)
		    day.append(int(cols[0]))
		    mon.append(cols[1])
		    yr.append(int(cols[2]))
		    mjd.append(float(cols[3]))
		    tJPL.append(float(cols[4]))
		    tacc.append(float(cols[5]))
		    f0.append(float(cols[6]))
		    f1.append(float(cols[8])*1.0E-15)
		    dm.append(float(cols[10]))
		    dm1.append(float(cols[11]))

		self.jtab = Table([day,mon,yr,mjd,tJPL,tacc,f0,f1,dm,dm1], 
		        names=('day', 'month', 'year', 'mjd', 'tJPL', 'tacc', 'f0', 'f1', 'dm', 'dm1'),
		        dtype=('i', 'S3', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

	def write_Crab_Tempo2_ephemeris(self,MJD,flag_write2file=True):
		diff = np.abs(self.jtab['mjd']-MJD)
		idx = diff.argmin()

		if diff[idx] > 18:
		    log.error('Too large a difference to closest ephemeris entry!')
		    sys.exit(1)
		    
		f0 = self.jtab[idx]['f0']
		f1 = self.jtab[idx]['f1']
		p0 = 1.0/f0
		p1 = -f1/(f0*f0)
		# This calculation of F2 comes from the C code at the end of the
		# explanatory notes for the Jodrell ephemeris
		f2 = 2.0*p1*p1/(p0*p0*p0)

		skel = """PSR B0531+21
RAJ 05:34:31.97232
DECJ +22:00:52.069
EPHEM DE200
F0 {0:.15f}
F1 {1:.15g}
C F2 computed as F2=2.0*P1*P1/(P0*P0*P0) according to JBO ephemeris
F2 {2:.15g}
PEPOCH {3:.15f}
TZRSITE @
TZRFRQ 0.0
TZRMJD {4:.15f}
START {5:.4f}
FINISH {6:.4f}
UNITS TDB
CLOCK TT(TAI)
DM {7:.5f}
DM1 {8:.5g}
DMEPOCH {9:.15f}
TRES {10:.3f}
"""

		# Note 5: The observed barycentric frequency and its first derivative
		# are quoted at the arrival time, using the DE200 ephemeris
		# So we use this value for both the period epoch and for TZRMJD
		# This is in TDB units
		mjd = self.jtab[idx]['mjd']
		pepoch = mjd + self.jtab[idx]['tJPL']/86400.0

		print(skel.format(f0,f1,f2,pepoch,pepoch,mjd-16,mjd+16,
		    self.jtab[idx]['dm'],self.jtab[idx]['dm1'],pepoch,self.jtab[idx]['tacc']))

		if flag_write2file:
			fname = 'crab_JBO_ephemeris_MJD%d.par' % MJD
			cmd = 'rm -f %s' % fname
			print(cmd);os.system(cmd)
			
			f = open(fname,'w')
			f.write(skel.format(f0,f1,f2,pepoch,pepoch,mjd-16,mjd+16,
			    self.jtab[idx]['dm'],self.jtab[idx]['dm1'],pepoch,self.jtab[idx]['tacc']))
			f.close()

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='make_crab_ephemeris_from_JBO.py',
		usage='make_crab_ephemeris_from_JBO.py MJD',
		description=code_description,
		epilog='',
		add_help=True,
		)
	parser.add_argument(
		'MJD',metavar='MJD',type=float,        
		help='MJD to extract ephemeris for')
	args = parser.parse_args()	
	
	ephem = Ephemeris_CrabPulsar()
	ephem.write_Crab_Tempo2_ephemeris(args.MJD)
