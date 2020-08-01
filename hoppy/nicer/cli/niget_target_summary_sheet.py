#!/usr/bin/env python

import os
import sys
import time
import datetime
import requests 
import argparse
from argparse import ArgumentParser
from bs4 import BeautifulSoup
import pandas as pd

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-01 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('niget_target_summary_sheet.py',
		usage='%(prog)s summary/segment [-u username] [-p password] [-o outcsvfile] ',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
(example) %(prog)s                [for "target segment sheet"] (default)
(example) %(prog)s --flag_summary [for "target summary sheet"] 

This script download the NICER "target summary sheet" or "target segment sheet"
as a csv format file from the NICER internal team webpage. 

The former is a summary of the observed target list, and the latter includes 
individual ObsIDs with corresponding observation dates, with the default output 
file names of nicer_target_segment_table.csv and nicer_target_summary_table.csv, 
respectively. The username and password are assumed to be set as environmental 
parameters. The following environmental parameters should be defined in your 
system.
	NICER_VIS_TEAM_WEBSITE : the webpage URL for the "target summary sheet"
	NICER_SEG_TEAM_WEBSITE : the webpage URL for the "target segment sheet"
	NICERPAGE_USERNAME : used for accessing the website 
	NICERPAGE_PASSWORD : used for accessing the website 
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--flag_summary', '-s', action='store_true',
		help='when flag_summary=True, the target_summary_sheet is downloaded instead of target_segement_details')
	parser.add_argument('--username', '-u', type=str, default=None, 
		help='username to access the NICER internal webpage')	
	parser.add_argument('--password', '-p', type=str, default=None, 
		help='password  to access the NICER internal webpage')
	parser.add_argument('--outcsvfile', '-o', type=str, default=None, 
		help='output csv filename')						
	return parser

def download_nicer_target_summary_sheet(username,password,outcsvfile):
	print("... Downloading of the NICER target summary sheet to a csv file.")

	if outcsvfile == None:
		outcsvfile = 'nicer_target_summary_table.csv'
	print('outcsvfile: {}'.format(outcsvfile))			

	if os.path.exists(outcsvfile):
		print("The output file has already existed, and is overwritten.")
		cmd = 'rm -f %s' % outcsvfile
		print(cmd);os.system(cmd)

	if os.getenv('NICER_VIS_TEAM_WEBSITE') is None:
		print("entironmental parameter NICER_VIS_TEAM_WEBSITE should be set.")
		quit()

	print("accessing to %s" % os.getenv('NICER_VIS_TEAM_WEBSITE'))
	req = requests.get(os.getenv('NICER_VIS_TEAM_WEBSITE'),auth=(username,password))
	if req.status_code != 200:
		sys.stderr.write('Error: problem on accessing the webpage.')
		quit()
	soup = BeautifulSoup(req.text,features="lxml")		
	tabs = soup.find_all('table')[4]	
	df = pd.read_html(str(tabs))[0]
	with open(outcsvfile, 'w') as f:
		f.write('# downlaoded %s \n' %  datetime.datetime.now().strftime('%Y%m%d_%H%M')) 
	df.to_csv(outcsvfile,mode='a')
	return outcsvfile

def download_nicer_target_segment_sheet(username,password,outcsvfile):
	print("... Downloading of the NICER target segment Details sheet to a csv file.")

	if outcsvfile == None:
		outcsvfile = 'nicer_target_segment_table.csv' 
	print('outcsvfile: {}'.format(outcsvfile))			

	if os.path.exists(outcsvfile):
		print("The output file has already existed, and is overwritten.")
		cmd = 'rm -f %s' % outcsvfile
		print(cmd);os.system(cmd)
		
	if os.getenv('NICER_SEG_TEAM_WEBSITE') is None:
		print("entironmental parameter NICER_SEG_TEAM_WEBSITE should be set.")
		quit()

	print("accessing to %s" % os.getenv('NICER_SEG_TEAM_WEBSITE'))
	req = requests.get(os.getenv('NICER_SEG_TEAM_WEBSITE'),auth=(username,password))
	if req.status_code != 200:
		sys.stderr.write('Error: problem on accessing the webpage.')
		quit()
	soup = BeautifulSoup(req.text,features="lxml")		
	tabs = soup.find_all('table')[1]
	df = pd.read_html(str(tabs))[0]
	df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
	with open(outcsvfile, 'w') as f:
		f.write('# downlaoded %s \n' %  datetime.datetime.now().strftime('%Y%m%d_%H%M')) 
	df.to_csv(outcsvfile,mode='a')
	return outcsvfile

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	if args.username is None:
		args.username = os.getenv("NICERPAGE_USERNAME")
	if args.password is None:
		args.password = os.getenv("NICERPAGE_PASSWORD")

	if args.username == None or args.password == None:
		sys.stderr.write('Error: username and password for the NICER webpage.\n')
		quit()
	print('username: {}'.format(args.username))
	print('password: {}'.format(args.password))	
	print('flag_summary: {}'.format(args.flag_summary))
	print('... False (default) for the segement details, True for the target summary')

	start_time = time.time()
	if args.flag_summary:
		outcsvfile = download_nicer_target_summary_sheet(
			username=args.username,
			password=args.password,
			outcsvfile=args.outcsvfile)
	else:
		outcsvfile = download_nicer_target_segment_sheet(
			username=args.username,
			password=args.password,
			outcsvfile=args.outcsvfile)
	elapsed_time = time.time() - start_time

if __name__=="__main__":
	main()