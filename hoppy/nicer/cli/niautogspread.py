#!/usr/bin/env python

import os 
import gspread
from oauth2client.service_account import ServiceAccountCredentials

scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']

credentials = ServiceAccountCredentials.from_json_keyfile_name(os.getenv('NICER_AUTO_GSPREAD_JSON'), scope)
gc = gspread.authorize(credentials)

webpage_address = "https://docs.google.com/spreadsheets/d/%s/edit#gid=0" % os.getenv('NICER_AUTO_GSPREAD_KEY')
print(webpage_address)

workbook = gc.open_by_key(os.getenv('NICER_AUTO_GSPREAD_KEY'))
print(workbook.title)

wks = workbook.sheet1
wks.update_acell('A1', 'Nyansuke!')
print(wks.acell('A1'))

