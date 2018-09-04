#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys 
import argparse

import math 

class SS433():
	def __init__(self,
		beta=0.2647,
		theta=20.92, # degree
		i=78.05, # degree
		P=162.375, # days
		t0=3563.23 # TJD
		):
		"""
		Eikenberry et al. 2001 ApJ, Table 1 
		http://adsabs.harvard.edu/abs/2001ApJ...561.1027E
		"""
		self.beta = beta
		self.theta = theta
		self.i = i
		self.P = P
		self.t0 = t0
		self.t0_MJD = self.t0 + 40000
		self.gamma = 1/math.sqrt(1-self.beta**2)
		print("beta: %.6f" % self.beta)
		print("theta: %.6f (degree) [%.6f]" % (self.theta,math.radians(self.theta)))
		print("i: %.6f (degree) [%.6f]" % (self.i,math.radians(self.i)))
		print("P: %.6f (days)" % self.P)
		print("t0: %.6f (TJD)" % self.t0)
		print("t0_MJD: %.6f (MJD)" % self.t0_MJD)		
		print("gamma: %.6f" % self.gamma)

	def get_redshift(self,MJD)	:
		phi = (MJD - self.t0_MJD)/self.P
		term1 = self.beta * math.sin(math.radians(self.theta)) * math.sin(math.radians(self.i)) * math.cos(phi)
		term2 = self.beta * math.cos(math.radians(self.theta)) * math.cos(math.radians(self.i))
		z1 = 1 - self.gamma * (1 + term1 + term2)
		z2 = 1 - self.gamma * (1 - term1 - term2)
		print(z1,z2,z1-z2)
		return z1, z2

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog=__file__,
		usage='%s MJD' % __file__,
		description='Calculate the Doppler shift of the SS433 jet precession.',
		epilog='',
		add_help=True,
		)
	parser.add_argument(
		'MJD',metavar='MJD',type=float,
		help='MJD.') 
	args = parser.parse_args()	

	ss433= SS433()
	redshift = ss433.get_redshift(args.MJD)
	print("input MJD: %.6f" % args.MJD)
	#print("redshift: %.6f" % redshift)
