#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 Sep. 7'
__version__ = '1.00'
"""
This script writes Tempo2-format PAR file from 
the Jodrell Bank Observatory (JBO) Crab monthly ephemeris parameters. 
Original code from Paul Ray. 
"""

class Ephemeris():
	def __init__(self):

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='make_crab_ephemeris_from_JBO.py',
		usage='make_crab_ephemeris_from_JBO.py MJD',
		description='A script to write Tempo2-format PAR file of the Crab pulsar from the JBP monthly ephemeris parameters.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'mjd',metavar='mjd',type=str,        
		help='Input MJD')
	args = parser.parse_args()	
	
