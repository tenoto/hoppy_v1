#!/usr/bin/env python

import os
import time
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
		description='(example) %(prog)s.py -o 1012010136 -y 2018_01 '
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--obsid', '-o', type=str, required=True, 
		help='input obsid (required)')	
	parser.add_argument('--yyyy_mm', '-y', type=str, required=True, 
		help='input yyyy_mm (required)')		
	return parser

def run_wget(obsid, yyyy_mm):
	cmd  = ''
	for subdir in ['auxil','xti','log']:
		cmd += 'wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np '
		cmd += '-R \'index*\' -erobots=off --retr-symlinks '
		cmd += 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/%s//%s/%s/;\n' % (yyyy_mm, obsid, subdir)
	print(cmd);os.system(cmd)

	return(cmd)



def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	print("--niwget--")
	print(" obsid: %s" % args.obsid)
	print(" yyyy_mm: %s" % args.yyyy_mm)

	start_time = time.time()
	cmd = run_wget(args.obsid, args.yyyy_mm)
	elapsed_time = time.time() - start_time

	if not os.path.exists(args.obsid):
		print("Error: file can not be download (check whether the yyyy_mm parameter is correct.)")
		quit()

	flog = '%s/wget.log' % args.obsid
	f = open(flog,'w')
	f.write('--niwget--\n')
	f.write(" yyyy_mm: %s\n" % args.yyyy_mm)	
	f.write(" obsid: %s\n" % args.obsid)	
	f.write(' requried time for download (sec): %.1f\n' % elapsed_time)
	f.write(' downloaded date: %s\n' % datetime.datetime.now())	
	f.write('\n')
	f.write(cmd)
	f.close()

if __name__=="__main__":
	main()