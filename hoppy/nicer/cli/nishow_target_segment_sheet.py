#!/usr/bin/env python

import os
import glob
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
	parser.add_argument('-s', action='store_true')		
	return parser

def show_nicer_target_segment_sheet(csvfile,target_name=None,obsid=None,flag_short=False):
	print("--input csvfile: {}".format(csvfile))
	print("--target name: {}".format(target_name))	
	print("--obsid: {}".format(obsid))		
	if target_name != None and obsid != None:
		print("Error: Please specify either of target_name or obsid options.")
		quit()

	df = pd.read_csv(csvfile,comment='#')
	df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
	if target_name != None:
		#print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']][df['Target Name'] == target_name])		
		for index, row in df.iterrows():
			flag_headas_repository = False
			flag_nicerteam_repository = False
			if target_name in row['Target Name']:
				file_path_headas_repository = '%s/nicer/data/obs/*/%s' % (os.getenv('HEADAS_REPOSITORY'),row['Observation ID'])
				if len(glob.glob(file_path_headas_repository)) > 0:
					flag_headas_repository = True
				file_path_nicerteam_repository = '%s/nicer/data/obs/*/%s' % (os.getenv('NICERTEAM_REPOSITORY'),row['Observation ID'])				
				if len(glob.glob(file_path_nicerteam_repository)) > 0:				
					flag_nicerteam_repository = True
				if flag_short:
					print(row['Observation ID'])
				else:
					print(row['Target Name'],row['Observation ID'],row['Start TimeUTC'],row['Good Expo[s]'],
						'HEADAS:%s' % flag_headas_repository,
						'NICERTEAM:%s' % flag_nicerteam_repository)
	elif obsid != None:		
		flag_headas_repository = False
		flag_nicerteam_repository = False
		file_path_headas_repository = '%s/nicer/data/obs/*/%s' % (os.getenv('HEADAS_REPOSITORY'),obsid)
		if len(glob.glob(file_path_headas_repository)) > 0:
			flag_headas_repository = True
		file_path_nicerteam_repository = '%s/nicer/data/obs/*/%s' % (os.getenv('NICERTEAM_REPOSITORY'),obsid)				
		if len(glob.glob(file_path_nicerteam_repository)) > 0:				
			flag_nicerteam_repository = True
		print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']][df['Observation ID'] == obsid],
			'HEADAS:%s' % flag_headas_repository,
			'NICERTEAM:%s' % flag_nicerteam_repository)		
	else:
		print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']])		

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	show_nicer_target_segment_sheet(args.csvfile,
		target_name=args.target_name,
		obsid=args.obsid,
		flag_short=args.s)

if __name__=="__main__":
	main()