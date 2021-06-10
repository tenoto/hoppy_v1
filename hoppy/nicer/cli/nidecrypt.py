#!/usr/bin/env python

import os
import argparse
from argparse import ArgumentParser

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-01 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('nidecrypt.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
		NICER data decrption 
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
	parser.add_argument('--obsid', '-o', type=str, default=None,
		help='target ObsID path')			
	parser.add_argument('--passphrase', '-p', type=str, default=None,
		help='decrypt passphrase')			
	return parser

def decrypt(obsid,passphrase):
	print("--obsid: {}".format(obsid))		

	cmd = "find %s -name '*.gpg' -print0 | xargs -n1 -0 gpg  --batch --yes --passphrase '%s';" % (obsid,passphrase)
	cmd += 'find %s -name "*.gpg" | xargs rm ' % obsid
	print(cmd);os.system(cmd)	

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	decrypt(args.obsid,args.passphrase)

if __name__=="__main__":
	main()

