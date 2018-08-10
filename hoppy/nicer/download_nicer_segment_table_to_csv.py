#!/usr/bin/env python

__name__    = 'download_nicer_segment_table_to_csv'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 April 5'

NICERPAGE = 'https://heasarc.gsfc.nasa.gov/docs/nicer/team_schedule/nicer_seg_team.html'

from optparse import OptionParser
import sys 
import requests
import pandas as pd
from bs4 import BeautifulSoup
import datetime 

parser = OptionParser()
parser.add_option("-u","--username",dest="username",
        action="store",help="NICER team username to access the webpage.",type="string")
parser.add_option("-p","--password",dest="password",
        action="store",help="NICER team password to access the webpage.",type="string")
parser.add_option("-o","--outcsv",dest="outcsvfile",default=None,
        action="store",help="Output CSV file",type="string")
(options, args) = parser.parse_args()

if options.outcsvfile == None:
	date = datetime.datetime.now()
	datestr = date.strftime('%Y%m%d_%H%M')
	options.outcsvfile = 'nicer_segment_table_v%s.csv' % datestr

if options.username == None or options.password == None:
	sys.stderr.write('Error: need username and password to access the NICER webpage.\n')
	quit()
req = requests.get(NICERPAGE, auth=(options.username,options.password))
if req.status_code != 200:
	sys.stderr.write('Error: problem on accessing the webpage.')
	quit()
soup = BeautifulSoup(req.text,'lxml')
tabs = soup.find_all('table')[1]
df = pd.read_html(str(tabs))[0]
df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)
df.to_csv(options.outcsvfile)



