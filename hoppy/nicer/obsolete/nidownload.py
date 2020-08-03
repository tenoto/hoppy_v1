__author__ = 'Teru Enoto'
__version__ = '1.03'
# v1.03 : modified for "click" environment 
# v1.02 : modified to use "download_nicer_signleObsID_data.py"
# v1.01 : original version

import os 
import sys 
import datetime 
import requests 
from bs4 import BeautifulSoup
import pandas as pd
import urllib.request

def url_is_alive(url):
    """
    Checks that a given URL is reachable.
    :param url: A URL
    :rtype: bool
    """
    request = urllib.request.Request(url)
    request.get_method = lambda: 'HEAD'

    try:
        urllib.request.urlopen(request)
        return True
    except urllib.request.HTTPError:
        return False

def wget_gsheet_nicer_download_target():
	# https://docs.google.com/spreadsheets/d/1Kj8lYz-ZvRYdKNRWOiDBeDk5wdNqRL6O7o4jTe9wKo4/edit#gid=645136886
	sheet_key = "1Kj8lYz-ZvRYdKNRWOiDBeDk5wdNqRL6O7o4jTe9wKo4"
	sheet_gid = "645136886"
	outfile = os.getenv("NICER_DOWNLOAD_SOURCE_LIST")
	cmd = 'wget "https://docs.google.com/spreadsheets/d/%s/export?format=csv&gid=%s" -O %s' % (sheet_key,sheet_gid,outfile)
	print(cmd);os.system(cmd)  
	with open(outfile, "r") as f:
		s = f.read()
	s = s.replace(",", " ")
	with open(outfile, "w") as f:
		f.write(s)

class NicerDownloader():
	def __init__(self):
		print("=== Object created: NicerDownloader ===")

	def get_target_segment_summary(self,username,password,outcsvfile):
		"""
		download NICER target segment summary table as csv file.
		"""
		print('--- method: %s ---' % sys._getframe().f_code.co_name)
		if username == None or password == None:
			sys.stderr.write('Error: need username and password to access the NICER webpage.\n')
			quit()
		print('username: {}'.format(username))
		print('password: {}'.format(password))	
		if outcsvfile == None:
			outcsvfile = 'nicer_segment_table_v%s.csv' % datetime.datetime.now().strftime('%Y%m%d_%H%M')
		print('outcsvfile: {}'.format(outcsvfile))			

		print('...accessing to {}'.format(os.getenv('NICER_SEG_TEAM_WEBSITE')))
		req = requests.get(os.getenv('NICER_SEG_TEAM_WEBSITE'),
			auth=(username,password))
		if req.status_code != 200:
			sys.stderr.write('Error: problem on accessing the webpage.')
			quit()
		soup = BeautifulSoup(req.text,'lxml')
		tabs = soup.find_all('table')[1]
		df = pd.read_html(str(tabs))[0]
		print(df)
		df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
		df.to_csv(outcsvfile)
		print("Done!")
		return outcsvfile

	def get_target_summary(self,username,password,outcsvfile):
		"""
		download NICER target summary table as csv file.
		"""
		print('--- method: %s ---' % sys._getframe().f_code.co_name)
		if username == None or password == None:
			sys.stderr.write('Error: need username and password to access the NICER webpage.\n')
			quit()
		print('username: {}'.format(username))
		print('password: {}'.format(password))	
		if outcsvfile == None:
			outcsvfile = 'nicer_target_table_v%s.csv' % datetime.datetime.now().strftime('%Y%m%d_%H%M')
		print('outcsvfile: {}'.format(outcsvfile))			

		print('...accessing to {}'.format(os.getenv('NICER_VIS_TEAM_WEBSITE')))
		req = requests.get(os.getenv('NICER_VIS_TEAM_WEBSITE'),
			auth=(username,password))
		if req.status_code != 200:
			sys.stderr.write('Error: problem on accessing the webpage.')
			quit()
		soup = BeautifulSoup(req.text,'lxml')
		print(len(soup.find_all('table')))
		tabs = soup.find_all('table')[4]
		df = pd.read_html(str(tabs))[0]
		print(df)
		df.to_csv(outcsvfile)
		print("Done!")
		return outcsvfile

	def get_single_obsid_data(self,obsid,yyyy_mm,decrypt):
		"""
		download a NICER single ObsID data set (public-->team)
		"""
		print('--- method: %s ---' % sys._getframe().f_code.co_name)		
		niobsid = NicerHeasarcObsID(obsid,yyyy_mm,gpg_decryption_password=decrypt)
		niobsid.wget_download()

	def show_target_segment_summary(self,csvfile,target_name=None):
		"""
		show a NICER target segment summary table in text format.
		"""
		print('--- method: %s ---' % sys._getframe().f_code.co_name)			
		df = pd.read_csv(csvfile)
		df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
		if target_name != None:
			print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']][df['Target Name'] == target_name])
		else:
			print(df[['Target Name','Observation ID','Start TimeUTC','Good Expo[s]']])

	def prepare_directories(self,dirpath_local_nicer_data='FTP/nicer/data/obs'):
		"""
		prepare directory setups.
		"""
		print('--- method: %s ---' % sys._getframe().f_code.co_name)				
		self.dirpath_local_nicer_data = dirpath_local_nicer_data
		if not os.path.exists(self.dirpath_local_nicer_data):
			cmd = 'mkdir -p %s' % self.dirpath_local_nicer_data
			print(cmd);os.system(cmd)
		if not os.path.exists('data'):
			cmd = 'ln -s FTP data'
			print(cmd);os.system(cmd)

	def set_download_source_list(self):
		"""
		set download source list from NICER_DOWNLOAD_SOURCE_LIST file.
		"""
		print('--- method: %s ---' % sys._getframe().f_code.co_name)			
		self.download_source_list = []
		for line in open(os.getenv('NICER_DOWNLOAD_SOURCE_LIST')):
			cols = line.split()
			if len(cols) == 0 or cols[0] == '#':
				continue
			self.download_source_list.append(cols[0])

	def get_multiple_obsid_data(self,username,password,decrypt):
		"""
		download multiple NICER ObsID datasets (public-->team).
		"""
		print('--- method: %s ---' % sys._getframe().f_code.co_name)			

		self.prepare_directories()
		self.set_download_source_list()
		print(self.download_source_list)	

		outcsvfile = self.get_target_segment_summary(username,password,outcsvfile=None)
		df = pd.read_csv(outcsvfile)
		df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)		

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

			#dirpath_yyyymm = '%s/FTP/nicer/data/obs/%s' % (os.getenv('HEASARC_REPOSITORY'), yyyy_mm)
			dirpath_yyyymm = './FTP/nicer/data/obs/%s' % yyyy_mm
			dirpath_obsid = '%s/%s' % (dirpath_yyyymm,obsid)
			if os.path.exists(dirpath_obsid):
				sys.stderr.write('target obsid %s is already existed. skipped.\n' % obsid)
				continue
			if not os.path.exists(dirpath_yyyymm):
				cmd = 'mkdir -p %s' % dirpath_yyyymm
				print(cmd);os.system(cmd)

			now = datetime.datetime.now().isoformat()
			message = '%s: %s %s %s\n' % (now,row['Target Name'],row['Observation ID'],row['Start TimeUTC'])
			sys.stdout.write(message)

			try:
				#cmd = 'download_nicer_singleObsID_data.py %s %s --gpg_decryption_password %s;' % (obsid,yyyy_mm,gpg_decryption_password)
				self.get_single_obsid_data(obsid,yyyy_mm,decrypt)				
				cmd = 'mv %s %s;\n' % (obsid, dirpath_yyyymm)
				print(cmd);os.system(cmd)
				message = '...downloaded %s\n' % obsid				
			except:
				message = '...error: download failed %s\n' % obsid
				sys.stderr.write(message)
				continue 
			sys.stdout.write(message)		
				

class NicerHeasarcObsID():
	def __init__(self,obsid,yyyy_mm,gpg_decryption_password=None):
		print("======== Object created: NicerHeasarcObsID =======")		
		self.obsid = obsid		
		self.yyyy_mm = yyyy_mm
		self.gpg_decryption_password = gpg_decryption_password
		print("obsid: {}".format(self.obsid))
		print("yyyy_mm: {}".format(self.yyyy_mm))
		print("gpg_decryption_password: {}".format(self.gpg_decryption_password))
		print("==================================")		

		if os.path.exists(self.obsid):
			sys.stderr.write('directory %s has already existed. Skipped.\n' % self.obsid)
			quit()

	def wget_download(self):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)				
		self.heasarc_public_path = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/%s/%s' % (self.yyyy_mm,self.obsid)
		self.heasarc_team_path   = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/.obs2validate/%s/%s' % (self.yyyy_mm,self.obsid)		
		self.attfile = 'ni%s.att.gz' % self.obsid
		if url_is_alive('%s/auxil/%s' % (self.heasarc_public_path,self.attfile)):
			self.flag_public = True
			print('target is the public website.')
			self.run_wget(self.heasarc_public_path)
			cmd = 'echo "%s: downloaded from the HEASARC public website." > %s/process_history.txt' % (
				datetime.datetime.now().strftime('%Y-%m-%d_%H:%M'),self.obsid)
			print(cmd);os.system(cmd)
		elif url_is_alive('%s/auxil/%s.gpg' % (self.heasarc_team_path,self.attfile)):
			self.flag_public = False 
			print('target is the team website.')			
			self.run_wget(self.heasarc_team_path)		
			if self.gpg_decryption_password != None:
				self.run_decrypt_gpt()	
			cmd = 'echo "%s: downloaded from the NICER team website." > %s/process_history.txt' % (
				datetime.datetime.now().strftime('%Y-%m-%d_%H:%M'),self.obsid)						
			print(cmd);os.system(cmd)				
		else:
			sys.stderr.write('[Error] No corresponding directory.\n')
			quit()
		print("Done!")

	def run_wget(self,target_path):
		subdir = 'xti'
		cmd  = 'wget -q -nH --no-check-certificate --cut-dirs=5 '
		cmd += '-r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks '
		cmd += '%s/%s;' % (target_path,subdir)		
		print(cmd);os.system(cmd)

	def run_decrypt_gpt(self):
		cmd  = 'decrypt_data.pl -d %s -p %s;' % (self.obsid,self.gpg_decryption_password)
		cmd += 'find %s -name "*.gpg" | xargs rm ' % self.obsid	
		print(cmd);os.system(cmd)













