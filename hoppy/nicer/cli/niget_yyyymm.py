#!/usr/bin/env python

import os
import argparse
from argparse import ArgumentParser
import pandas as pd

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-01 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('niget_yyyymm.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
		help description
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--csvfile', '-c', type=str, default="nicer_target_segment_table.csv", 
		help='csvfile')	
	parser.add_argument('--obsid', '-o', type=str, default=None,
		help='target ObsID (default=None)')			
	return parser

def get_yyyymm(csvfile,obsid):
	print("--input csvfile: {}".format(csvfile))
	print("--obsid: {}".format(obsid))		

	df = pd.read_csv(csvfile,comment='#')
	df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
	str_start_timeutc = df['Start TimeUTC'][df['Observation ID'] == obsid]
	print(str_start_timeutc)

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	get_yyyymm(args.csvfile,obsid=args.obsid)

if __name__=="__main__":
	main()