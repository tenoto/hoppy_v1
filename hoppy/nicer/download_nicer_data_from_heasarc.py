#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 August 10'
__version__ = '1.01'

import os 
import sys 
import argparse 
import datetime
import pandas as pd 

class DownloadManager():
	def __init__(self,dir_input='script/input'):
		self.dir_input = dir_input

		if os.getenv('HEASARC_REPOSITORY') == "":
			sys.stderr.write('HEASARC_REPOSITORY is not set.\n')
			quit()
		os.chdir(os.getenv('HEASARC_REPOSITORY'))

	def download_observation_segment(self,username,password,flag_download_newtable=True):
		self.username = username
		self.password = password
		self.flag_download_newtable = flag_download_newtable		

		if not os.path.exists(self.dir_input):
			cmd = 'mkdir -p %s' % (self.dir_input)
			print(cmd);os.system(cmd)

		if self.flag_download_newtable:
			date = datetime.datetime.now()
			datestr = date.strftime('%Y%m%d_%H%M')
			self.fname_segment_table = '%s/nicer_segment_table_v%s.csv' % (self.dir_input,datestr)

			print("...downloading the latest segment table...")
			cmd = 'rm -f %s/nicer_segment_table_v*.csv\n' % self.dir_input
			cmd += 'download_nicer_segment_table_to_csv.py '
			cmd += '--username %s ' % self.username
			cmd += '--password %s ' % self.password
			cmd += '--outcsv %s ' % self.fname_segment_table
			print(cmd);os.system(cmd)
		else:
			csvlst = glob.glob('%s/nicer_segment_table_v*.csv' % self.dir_input)
			if len(csvlst) == 0:
				sys.stderr.write('no nicer_segment_table_v*.csv\n')
				quit()
			else:
				self.fname_segment_table = csvlst[0]
		sys.stdout.write('fname_segment_table:%s\n' % self.fname_segment_table)

	def prepare_download_source_list(self):
		cmd  = 'rm -f %s/download_source_list.txt;\n' % self.dir_input
		cmd += 'cp %s %s' % (os.getenv('NICER_SOURCE_LIST_FOR_DOWNLOAD'),self.dir_input)
		print(cmd);os.system(cmd)

		self.download_source_list = []
		for line in open(os.getenv('NICER_SOURCE_LIST_FOR_DOWNLOAD')):
			cols = line.split()
			if len(cols) == 0 or cols[0] == '#':
				continue
			self.download_source_list.append(cols[0])
		print(self.download_source_list)	

	def prepare_directories(self,dirpath_local_nicer_data='FTP/nicer/data/obs'):
		self.dirpath_local_nicer_data = dirpath_local_nicer_data
		if not os.path.exists(self.dirpath_local_nicer_data):
			cmd = 'mkdir -p %s' % self.dirpath_local_nicer_data
			print(cmd);os.system(cmd)
		if not os.path.exists('data'):
			cmd = 'ln -s FTP data'
			print(cmd);os.system(cmd)

	def downalod_source(self):
		df = pd.DataFrame.from_csv(self.fname_segment_table)
		df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)

		self.fname_logfile = '%s/download.log' % self.dir_input
		f_logfile = open(self.fname_logfile,'a')
		for index, row in df.iterrows():
			if not row['Target Name'] in self.download_source_list:
				continue 
			if not row['Process State'] in ['DELV','DELVV']:
				continue 
			if row['Observation ID'][0] == '0':
				continue

			try:
				dtime = datetime.datetime.strptime(row['Start TimeUTC'], "%Y-%m-%dT%H:%M:%S")
			except:
				print("... skip (ObsID=%s) because Start TimeUTC is blank." % obsid)
				continue 
			yyyy_mm = "{}_{:02d}".format(dtime.year, dtime.month)
			obsid = str(row['Observation ID'])
			tarfile = '%s.tar' % obsid
			target_path = '%s/nicer/data/obs/%s/%s' % (os.getenv('HEASARC_REMOTE_FTPDIR'), yyyy_mm, obsid)
			dir_move_to = '%s/FTP/nicer/data/obs/%s' % (os.getenv('HEASARC_REPOSITORY'), yyyy_mm)

			dir_final_path = '%s/%s' % (dir_move_to,obsid)
			if os.path.exists(dir_final_path):
				sys.stderr.write('target obsid %s is already existed. skipped.\n' % obsid)
				continue

			if not os.path.exists(dir_move_to):
				cmd = 'mkdir -p %s' % dir_move_to
				print(cmd);os.system(cmd)

			now = datetime.datetime.now().isoformat()
			message = '%s: %s %s %s\n' % (now,row['Target Name'],row['Observation ID'],row['Start TimeUTC'])
			sys.stdout.write(message)
			f_logfile.write(message)
			try:
				#cmd = ''
				#for subdir in ['auxil','log','xti']:
				#	cmd += 'wget -q -nH --no-check-certificate --cut-dirs=5 '
				#	cmd += '-r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks '
				#	cmd += '%s/%s;' % (target_path,subdir)
				subdir = 'xti'
				cmd  = 'wget -q -nH --no-check-certificate --cut-dirs=5 '
				cmd += '-r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks '
				cmd += '%s/%s;' % (target_path,subdir)		
				cmd += 'mv %s %s;\n' % (obsid, dir_move_to)
				print(cmd);os.system(cmd)
				f_logfile.write(cmd)
				message = '...downloaded %s\n' % obsid
			except:
				message = '...error: download faield %s\n' % obsid
				sys.stderr.write(message)
			sys.stdout.write(message)
			f_logfile.write(message)
			#exit()
		f_logfile.close()


if __name__=="__main__":

	sys.stdout.write('... download a NICER observagtion from the heasarc archive ...\n')

	parser = argparse.ArgumentParser(
		prog='download_nicer_data_from_heasarc.py',
		usage='download_nicer_data_from_heasarc.py [username] [password] [newtable]',
		description='Download script of a NICER data from the heasarc archive.',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'-u','--username',metavar='username',type=str,     
		help='NICER team username to access the webpage.')
	parser.add_argument(
		'-p','--password',metavar='password',type=str,     
		help='NICER team password to access the webpage.')	
	parser.add_argument(
		'-n', '--newtable',action='store_true',dest='flag_download_newtable',
		default=True,
		help='flag to download a new segment table.')	
	args = parser.parse_args()	
	
	dm = DownloadManager()
	dm.download_observation_segment(args.username,args.password,args.flag_download_newtable)
	dm.prepare_download_source_list()
	dm.prepare_directories()
	dm.downalod_source()