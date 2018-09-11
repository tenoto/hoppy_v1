#!/usr/bin/env python

import os 
import sys 
import glob
import yaml 
import argparse

#import hoppy.nicer.nipipeline as nip

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='run_nixspec_multi_observations.py',
		usage='run_nixspec_multi_observations.py obsid_path_str [--outdir (default=nispec)] [--flag_evt]',
		description='run series of xspec fitting for NICER data.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'obsid_path_str',metavar='obsid_path_str',type=str,        
		help='Input ObsID directory for glob. e.g., "../../../pipe/v180809_ws1/niout/o1u20cut/10500501*"')
	parser.add_argument(
		'--outdir',metavar='--outdir',type=str,default="nispec",
		help='Output directory (default=nispec).')
	parser.add_argument('--flag_evt', action='store_true', 
		help='flag to copy of event files.') 	
	parser.add_argument(
		'--rmffile',metavar='--rmffile',type=str,
		default="%s/nicer_v1.02.rmf" % os.getenv('NICER_RESP_PATH'),
		help='rmf file.')
	parser.add_argument(
		'--arffile',metavar='--arffile',type=str,
		default="%s/ni_xrcall_onaxis_v1.02_wo_detid14_34.arf" % os.getenv('NICER_RESP_PATH'),
		help='arf file.')	
	parser.add_argument(
		'--modelxcm',metavar='--modelxcm',type=str,
		default="%s/run_nixspec_multi_observations_bb.xcm" % os.getenv('NICER_SOFT_PATH'),
		help='model xcm.')	
	parser.add_argument(
		'--inputyaml',metavar='--inputyaml',type=str,
		default="%s/run_nixspec_multi_observations.yaml" % os.getenv('NICER_SOFT_PATH'),
		help='input yaml file.')				
	args = parser.parse_args()	

	print(args)

	# preapre directory 
	dir_spec = '%s/pha' % args.outdir
	cmd  = 'prepare_nicer_xspec_directory.py '
	cmd += '"%s" ' % args.obsid_path_str
	cmd += '--outdir %s' % dir_spec 
	if args.flag_evt:
		cmd += '--flag_evt '
	print(cmd);os.system(cmd)

	# make csvfile 
	dir_param = '%s/input' % args.outdir
	cmd = 'mkdir -p %s' % dir_param
	print(cmd);os.system(cmd)
	fname_filelist2csv_txt = '%s/filelist2csv.txt' % dir_param
	f = open(fname_filelist2csv_txt,'w')
	for obsid_path in glob.glob('%s/*' % dir_spec):
		obsid = os.path.basename(obsid_path)
		src_pha_list = glob.glob('%s/%s/ni%s_*_clscr_gtisel.pha' % (dir_spec,obsid,obsid))
		if len(src_pha_list) == 0:
			continue
		src_pha = src_pha_list[0]
		bgd_pha_list = glob.glob('%s/%s/ni%s_*_BGMod_3C50.pha' % (dir_spec,obsid,obsid))
		if len(bgd_pha_list) == 0:
			continue 
		bgd_pha = bgd_pha_list[0]
		dump = '%s,%s,%s,%s,%s' % (
			obsid,
			src_pha,
			bgd_pha,
			args.rmffile,args.arffile)
		print(dump)
		f.write(dump+'\n')
	f.close()

	cmd = 'cp %s %s' % (args.modelxcm, dir_param)
	print(cmd);os.system(cmd)

	dir_fit = '%s/fit' % args.outdir
	fname_yaml_output = '%s/run_nixspec_multi_observations.yaml' % dir_param
	param = yaml.load(open(args.inputyaml))
	param['outdir'] = dir_fit
	param['modelxcm'] = '%s/%s' % (dir_param,os.path.basename(args.modelxcm))
	with open(fname_yaml_output, 'w') as outfile:		
		yaml.dump(param, outfile, default_flow_style=False)

	fname_fit_csv = '%s/xspec_multifit.csv' % dir_param
	cmd  = 'make_csv2xspec.py '
	cmd += '%s ' % fname_filelist2csv_txt
	cmd += '%s ' % fname_yaml_output
	cmd += '%s ' % fname_fit_csv
	print(cmd);os.system(cmd)

	dump  = "#!/bin/sh -f\n"
	dump += 'rm -rf %s\n' % dir_fit 
	dump += 'xspec_multi_observations.py '
	dump += '%s ' % fname_fit_csv
	dump += '%s ' % fname_yaml_output
	fname_run_sh = '%s/run_fit.sh' % dir_param
	f = open(fname_run_sh,'w')
	f.write(dump)
	f.close()
	cmd = 'chmod +x %s' % fname_run_sh
	print(cmd);os.system(cmd)
	print("==========================")
	print("run following command")
	print("==========================")	
	print("%s" % fname_run_sh)



