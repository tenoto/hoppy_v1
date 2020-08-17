#!/usr/bin/env python

import os
import sys
import time
import pandas as pd
import datetime
import argparse

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-01 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('niwget',
		usage='%(prog)s -o obsid -y yyyy_mm ',
		description="""
(example) %(prog)s.py -o 1012010136 -y 2018_01 
(example) %(prog)s.py -t "comet"
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--obsid', '-o', type=str, 
		help='input obsid)')	
	parser.add_argument('--yyyy_mm', '-y', type=str, 
		help='input yyyy_mm')	
	parser.add_argument('--target_name', '-t', type=str, default=None,
		help='target source name to be shown (default=None)')	
	parser.add_argument('--csvfile', '-c', type=str, default="nicer_target_segment_table.csv", 
		help='csvfile')		
	parser.add_argument('--heasarc_repository', '-r', type=str, default=os.getenv('HEADAS_REPOSITORY'),
		help='Heasarc repository directory. If this option is specified, the files are moved here.')	

	return parser

def run_wget(obsid, yyyy_mm, heasarc_repository=None):
	print("\n[niwget] %s" % (sys._getframe().f_code.co_name))

	if os.path.exists(obsid):
		print("Error: directory has already existed.")
		return -1

	if heasarc_repository != None:
		if os.path.exists('%s/nicer/data/obs/%s/%s' % (heasarc_repository,yyyy_mm,obsid)):
			print("Error: directory has already existed. skip.")
			return -1

	start_time = time.time()
	cmd  = ''
	for subdir in ['auxil','xti','log']:
		cmd += 'wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np '
		cmd += '-R \'index*\' -erobots=off --retr-symlinks '
		cmd += 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/%s//%s/%s/;\n' % (yyyy_mm, obsid, subdir)
	print(cmd);os.system(cmd)
	elapsed_time = time.time() - start_time

	if not os.path.exists(obsid):
		print("Error: file can not be download (check whether the yyyy_mm parameter is correct.)")
		return -1

	flog = '%s/wget.log' % obsid
	f = open(flog,'w')
	f.write('--niwget--\n')
	f.write(" yyyy_mm: %s\n" % yyyy_mm)	
	f.write(" obsid: %s\n" % obsid)	
	f.write(' requried time for download (sec): %.1f\n' % elapsed_time)
	f.write(' downloaded date: %s\n' % datetime.datetime.now())	
	f.write('\n')
	f.write(cmd)
	f.close()

	if heasarc_repository != None:
		outdir = '%s/nicer/data/obs/%s' % (heasarc_repository,yyyy_mm)
		if not os.path.exists(outdir):
			cmd = 'mkdir -p %s' % outdir
			print(cmd);os.system(cmd)
		if os.path.exists('%s/%s' % (outdir,obsid)):
			cmd = 'rm -rf %s/%s' % (outdir,obsid)
			print(cmd);os.system(cmd)
		cmd = 'mv %s %s' % (obsid,outdir)
		print(cmd);os.system(cmd)

	return(cmd)

def download_single_obsid(args):
	print("\n[niwget] %s" % (sys._getframe().f_code.co_name))

	print("--niwget--")
	print(" obsid: %s" % args.obsid)
	print(" yyyy_mm: %s" % args.yyyy_mm)
	run_wget(args.obsid, args.yyyy_mm, heasarc_repository=args.heasarc_repository)

def download_multiple_obsid(args):
	print("\n[niwget] %s" % (sys._getframe().f_code.co_name))

	print("target_name: %s" % args.target_name)
	if args.csvfile == None or not os.path.exists(args.csvfile):
		print("Error: pleaes add --csvfile for the nicer observation list.")
		exit()

	df = pd.read_csv(args.csvfile,comment='#')
	df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
	for index, row in df.iterrows():
		if args.target_name in row['Target Name']:
			print(row['Target Name'],row['Observation ID'],row['Start TimeUTC'],row['Good Expo[s]'])
			yyyy, mm, others = row['Start TimeUTC'].split('-')
			yyyy_mm = yyyy + "_" + mm
			run_wget(row['Observation ID'], yyyy_mm, heasarc_repository=args.heasarc_repository)

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	if args.target_name == None:
		download_single_obsid(args)
	else:
		download_multiple_obsid(args)

if __name__=="__main__":
	main()