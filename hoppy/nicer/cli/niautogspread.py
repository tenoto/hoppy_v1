#!/usr/bin/env python

import os 
import pandas as pd
import gspread
from gspread_dataframe import set_with_dataframe, get_as_dataframe
from oauth2client.service_account import ServiceAccountCredentials

scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']

credentials = ServiceAccountCredentials.from_json_keyfile_name(os.getenv('NICER_AUTO_GSPREAD_JSON'), scope)
gc = gspread.authorize(credentials)

webpage_address = "https://docs.google.com/spreadsheets/d/%s/edit#gid=0" % os.getenv('NICER_AUTO_GSPREAD_KEY')
#print(webpage_address)

workbook = gc.open_by_key(os.getenv('NICER_AUTO_GSPREAD_KEY'))
#print(workbook.title)

# ==========
dictkey = {
	'Observatory':'A',
	'ObsID':'B',
	'Note':'C',
	'Download':'D',
	'Pipeline':'E',
	'Spectrum':'D'}

worksheet = workbook.worksheet('process')

df = pd.DataFrame(worksheet.get_all_records())
#df = get_as_dataframe(worksheet)
for index, row in df.iterrows():
	df['Download'][index] = 'HEASARC'
	cellid = '%s%d' % (dictkey['Download'],index+2)
	worksheet.format(cellid, {"backgroundColor": {"red": 0.0,"green": 1.0, "blue": 0.0}})

#print(df.reset_index(drop=True))
#print(df.reset_index(drop=True, inplace=True))
set_with_dataframe(worksheet, df.reset_index(drop=True))


#	worksheet.update([df.columns.values.tolist()]+df.values.tolist())
#worksheet.update('A3', 'Bingo!')


