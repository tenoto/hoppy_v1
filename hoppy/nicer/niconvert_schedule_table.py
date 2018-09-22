#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 September 22'
__version__ = '1.00'

import os 
import sys 
import argparse
import hoppy.timing.xtime as xt

def convert_nicer_schedule_file(fname_schedule):
	if not os.path.exists(fname_schedule):
		sys.stderr.write(fname_schedule)
		exit()

	dump  = 'Target     Visibility Start      Visibility End       Span\n'       
	dump += 'Name       (UTC)                 (UTC)                (seconds)\n'  
	dump += '--------   -----------------     -----------------    -------\n'
	for line in open(fname_schedule):
		cols = line.split()
		if len(cols) == 0:
			continue
		if cols[0] == 'Target' or cols[0] == 'Name' or cols[0] == '------------------':
			continue
	 	target_name = cols[0]
		visibility_start = cols[2]
		visibility_end = cols[3]
		visibility_span = float(cols[4])
		#print(visibility_start,visibility_end,visibility_span)

		#visibility_start
		t0 = xt.XTime(
			visibility_start.replace('T',':').replace('-',':'),
			'yday','utc')
		t1 = xt.XTime(
			visibility_end.replace('T',':').replace('-',':'),
			'yday','utc')		
		dump += '%s     %s     %s  %.2f\n' % (target_name, t0.get_utc(format='isot'),t1.get_utc(format='isot'),visibility_span)
	print(dump)




if __name__=="__main__":

	parser = argparse.ArgumentParser(
		prog='niconvert_schedule_table.py',
		usage='niconvert_schedule_table.py schedule.file',
		description='A script to convert NICER schedule table to UTC ISOT format.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'fname_schedule',metavar='fname_schedule',type=str,        
		help='NICER schedule file')
	args = parser.parse_args()	
	print(args)

	convert_nicer_schedule_file(args.fname_schedule)
	
