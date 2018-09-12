#!/usr/bin/env python

import os 
import sys 
import glob
import yaml 
import argparse
import pandas as pd 
#import hoppy.nicer.nipipeline as nip

if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='plot_nixspec_multi_observations.py',
		usage='plot_nixspec_multi_observations.py fit.csv input.yaml',
		description='plot series of xspec fitting for NICER data.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'fitcsv',metavar='fitcsv',type=str,        
		help='fit csvfile')
	parser.add_argument(
		'inputyaml',metavar='inputyaml',type=str,
		help='input yaml file.')		
	parser.add_argument(
		'--rateband_nsel',metavar='--rateband_nsel',type=int,default=0,
		help='rateband selected number')
	args = parser.parse_args()	

	print(args)

	param = yaml.load(open(args.inputyaml))	
	df = pd.read_csv(args.fitcsv)	

	emin_rate, emax_rate = param['ratebands'][args.rateband_nsel]		
	ycol_name = 'rate_%sto%skeV' % (str(emin_rate).replace('.','p'),str(emax_rate).replace('.','p'))
	ycolerr_name = '%s_err' % ycol_name
	plot_xmin = int(min(df['MJD_DATEOBS']-20))
	plot_xmax = int(max(df['MJD_DATEOBS']+20))	
	plot_ymax = max(df[ycol_name])*1.5
	plot_ymin = -max(df[ycol_name])*0.02
	
	fname_fluxplot_yaml_input = '%s/run_nixspec_multi_observations_fluxplot.yaml' % os.getenv('NICER_SOFT_PATH')
	fname_fluxplot_yaml = '%s/run_nixspec_multi_observations_fluxplot.yaml' % os.path.dirname(args.inputyaml)
	param_plot = yaml.load(open(fname_fluxplot_yaml_input))
	param_plot['datagroup'][0]['file'] = '%s/xspec_multifit_fit.csv' % param['outdir']
	param_plot['datagroup'][0]['label'] = "NICER (%.1f-%.1f keV)" % (emin_rate,emax_rate)
	param_plot['datagroup'][0]['ycolumn'] = ycol_name
	param_plot['datagroup'][0]['ycolumn_error'] = ycolerr_name	
	#param_plot['ycolumn'] = ycol_name
	#param_plot['ycolumn_error'] = ycolerr_name
	param_plot['outpdf'] = '%s/run_nixspec_multi_observations_fluxplot.pdf' % os.path.dirname(args.fitcsv)
	param_plot['xlim'] = [plot_xmin,plot_xmax]
	param_plot['ylim'] = [plot_ymin,plot_ymax]
	with open(fname_fluxplot_yaml, 'w') as outfile:		
		yaml.dump(param_plot, outfile, default_flow_style=False)

	cmd  = 'plot_x1y1_multigroup.py '
	cmd += '%s ' % fname_fluxplot_yaml
	print(cmd);os.system(cmd)
