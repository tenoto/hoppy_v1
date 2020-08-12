#!/usr/bin/env python

import os
import yaml
import pandas as pd
import argparse

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2020-08-12 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('niauto.py',
		usage='%(prog)s -o obsid',
		description="""
(example) %(prog)s.py -o 1012010136
(example) %(prog)s.py -t "comet"
"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('obsid',metavar='obsid',type=str,help='input obsid')
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command')
#	parser.add_argument('--target_name', '-t', type=str, default=None,
#		help='target source name to be shown (default=None)')	
	parser.add_argument('--csvfile', '-c', type=str, default="nicer_target_segment_table.csv", 
		help='csvfile')		
	parser.add_argument('--setupfile', '-s', type=str, default="{}/hoppy/nicer/cli/setup_template.yaml".format(os.getenv('HOPPY_PATH')), 
		help='setup yaml file')		
	parser.add_argument('--heasarc_repository', '-r', type=str, default=os.getenv('HEADAS_REPOSITORY'),
		help='Heasarc repository directory. If this option is specified, the files are moved here.')	

	return parser

def niauto(args):
	if not os.path.exists(args.csvfile):
		cmd = 'niget_target_segment_sheet.py'
		print(cmd);os.system(cmd)			

	df = pd.read_csv(args.csvfile,comment='#')
	df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
	str_start_timeutc = df['Start TimeUTC'][df['Observation ID'] == args.obsid]
	yyyy = str(str_start_timeutc).split('-')[0].split()[-1]
	mm = str(str_start_timeutc).split('-')[1]	
	yyyy_mm = "%s_%s" % (yyyy,mm)

	datadir = '%s/nicer/data/obs/%s/%s' % (os.getenv('HEADAS_REPOSITORY'),yyyy_mm,args.obsid)
	if not os.path.exists(datadir):
		cmd = 'niwget.py -o %s -y %s -r %s' % (args.obsid,yyyy_mm,args.heasarc_repository)
		print(cmd);os.system(cmd)

	outdir = 'out/%s' % args.obsid
	if os.path.exists(outdir):
		print("output directory %s has already existed." % outdir)
		exit()

	cmd = 'mkdir -p %s;' % (outdir)
	print(cmd);os.system(cmd)

	param = yaml.load(open(args.setupfile))
	param['output_directory'] = outdir
	input_setup_yaml = '%s/input_setup.yaml' % outdir
	with open(input_setup_yaml, 'w') as file:
	    yaml.dump(param, file)
	input_obsid_lst = '%s/input_obsid.lst' % outdir
	f = open(input_obsid_lst,'w')
	f.write('%s\n' % args.obsid)
	f.close()

	cmd  = 'nipipeline.py '
	cmd += '--setup_yamlfile %s ' % input_setup_yaml
	cmd += '--obsid_lstfile %s ' % input_obsid_lst
	print(cmd);os.system(cmd)

	print("finished...\n")

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	niauto(args)

if __name__=="__main__":
	main()