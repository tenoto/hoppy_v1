#!/usr/bin/env python

import os
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
		description='(example) %(prog)s.py -o 1012010152 -y 2018_05 '
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

def main(args=None):
	"""
	Main entry point for your project.
	Args:
		args : list
		A of arguments as if they were input in the command line. 
		Leave it None to use sys.argv.
    """
	parser = get_parser()
	args = parser.parse_args(args)

	print("--niwget--")
	print(" obsid: %s" % args.obsid)
	print(" yyyy_mm: %s" % args.yyyy_mm)
	run_wget(args.obsid, args.yyyy_mm)

if __name__=="__main__":
	main()