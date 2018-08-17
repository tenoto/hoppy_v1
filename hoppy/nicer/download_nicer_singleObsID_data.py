#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__  = 'Teru Enoto'
__date__    = '2018 August 17'
__version__ = '1.00'
# v1.00 : original version 

import os 
import sys 
import argparse 
import urllib2

def exists(location):
	print("%s exists?" % location)
	request = urllib2.Request(location)
	request.get_method = lambda : 'HEAD'
	try:
		response = urllib2.urlopen(request)
		sys.stdout.write('Yes, file exists.\n')
		return True
	except urllib2.HTTPError:
		sys.stdout.write('No, file does not exist.\n')		
		return False

class HeasarcObsID():
	def __init__(self,obsid,yyyy_mm,gpg_decryption_password=None):
		self.obsid = obsid		
		self.yyyy_mm = yyyy_mm
		print("==================================")
		print("obsid: %s" % self.obsid)
		print("yyyy_mm: %s" % self.yyyy_mm)

		self.gpg_decryption_password = gpg_decryption_password
		print("gpg_decryption_password: %s" % self.gpg_decryption_password)
		print("==================================")		

		if os.path.exists(self.obsid):
			sys.stderr.write('directory %s has already existed. Skipped.\n' % self.obsid)
			quit()

	def wget_download(self):
		self.heasarc_public_path = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/%s/%s' % (self.yyyy_mm,self.obsid)
		self.heasarc_team_path   = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/.obs2validate/%s/%s' % (self.yyyy_mm,self.obsid)		
		self.attfile = 'ni%s.att.gz' % self.obsid
		if exists('%s/auxil/%s' % (self.heasarc_public_path,self.attfile)):
			self.flag_public = True
			sys.stdout.write('target is the public website.\n')			
			self.run_wget(self.heasarc_public_path)
			cmd = 'touch %s/downloaded_from_public_website.txt' % self.obsid
			print(cmd);os.system(cmd)
		elif exists('%s/auxil/%s.gpg' % (self.heasarc_team_path,self.attfile)):
			self.flag_public = False 
			sys.stdout.write('target is the team website.\n')			
			self.run_wget(self.heasarc_team_path)		
			if self.gpg_decryption_password != None:
				self.run_decrypt_gpt()	
			cmd = 'touch %s/downloaded_from_team_website.txt' % self.obsid
			print(cmd);os.system(cmd)				
		else:
			sys.stderr.write('[Error] No corresponding directory.\n')
			quit()

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

if __name__=="__main__":

	sys.stdout.write('... download a NICER single ObsID data set (public-->team)...\n')

	parser = argparse.ArgumentParser(
		prog='download_nicer_singleObsID_data.py',
		usage='download_nicer_singleObsID_data.py obsid yyyy_mm [--gpg_decryption_password]',
		description='Download a NICER single ObsID data set (public-->team)',
		epilog='',
		add_help=True,
		)

	parser.add_argument(
		'obsid',metavar='obsid',type=str,     
		help='ObsID')	
	parser.add_argument(
		'yyyy_mm',metavar='yyyy_mm',type=str,     
		help='yyyy_mm')	
	parser.add_argument(
		'-p','--gpg_decryption_password',metavar='password',type=str,     
		default=None,
		help='NICER team password to access the webpage.')	
	args = parser.parse_args()	

	obsid = HeasarcObsID(args.obsid,args.yyyy_mm,gpg_decryption_password=args.gpg_decryption_password)
	obsid.wget_download()
