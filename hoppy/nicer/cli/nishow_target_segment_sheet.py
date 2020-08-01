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
	parser = argparse.ArgumentParser('nishow_target_segment_sheet.py',
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
	parser.add_argument('--target_name', '-t', type=str, default=None,
		help='target source name to be shown (default=None)')		
	parser.add_argument('--obsid', '-o', type=str, default=None,
		help='target ObsID (default=None)')			
	return parser

def show_nicer_target_segment_sheet(csvfile,target_name=None,obsid=None):
	print("--input csvfile: {}".format(csvfile))
	print("--target name: {}".format(target_name))	
	print("--obsid: {}".format(obsid))		
	if target_name != None and obsid != None:
		print("Error: Please specify either of target_name or obsid options.")
		quit()

	df = pd.read_csv(csvfile,comment='#')
	df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
	if target_name != None:
		print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']][df['Target Name'] == target_name])
	elif obsid != None:		
		print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']][df['Observation ID'] == obsid])		
	else:
		print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']])		

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	show_nicer_target_segment_sheet(args.csvfile,
		target_name=args.target_name,
		obsid=args.obsid)

if __name__=="__main__":
	main()