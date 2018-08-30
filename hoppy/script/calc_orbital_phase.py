#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import sys 
import argparse

import math 
import sympy 

class BinaryOrbit():
	def __init__(self,a0,a1,a2=0.0):
		self.a0 = a0
		self.a1 = a1
		self.a2 = a2
		print("a0: %.6f" % self.a0)
		print("a1: %.6f" % self.a1)
		print("a2: %.6f" % self.a2)			

	def get_orbital_phase(self,MJD)	:
		n_param = sympy.Symbol('n_param')
		phase_solutions = sympy.solve(self.a0 + self.a1 * n_param + self.a2 * n_param ** 2 - MJD, n_param)
		phase_solutions_plus = []
		for i in phase_solutions:
			if i > 0.0:
				phase_solutions_plus += [i]
		return phase_solutions_plus[0]

	def get_mjd(self,n_cycle_phase):
		return self.a0 + self.a1 * n_cycle_phase + self.a2 * n_cycle_phase ** 2 


if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog=__file__,
		usage='%s MJD a0 a1 a2=0' % __file__,
		description='Calculate the orbital phase from an input MJD value with coefficient of an orbital ephemeris.',
		epilog='',
		add_help=True,
		)
	parser.add_argument(
		'MJD',metavar='MJD',type=float,
		help='MJD.') 
	parser.add_argument(
		'a0',metavar='a0',type=float,
		help='Sufficient of 0th order (constant) of the orbital ephemeris.') 	
	parser.add_argument(
		'a1',metavar='a1',type=float,
		help='Sufficient of 1st order (linear) of the orbital ephemeris.') 		
	parser.add_argument(
		'a2',metavar='a2',type=float,default=0.0,
		help='Sufficient of 1st order (linear) of the orbital ephemeris.') 			
	args = parser.parse_args()	

	orbit = BinaryOrbit(args.a0,args.a1,args.a2)
	orbital_cycle = orbit.get_orbital_phase(args.MJD)
	print("input MJD: %.6f" % args.MJD)
	print("orbital_cycle: %.6f" % orbital_cycle)
	print("orbital_phase: %.6f" % math.modf(orbital_cycle)[0])
	mjd_crosscheck = orbit.get_mjd(orbital_cycle)
	print("MJD (sanity check): %.6f" % mjd_crosscheck)
