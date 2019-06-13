#!/usr/bin/env python

import os 
import sys 
if len(sys.argv) != 3:
	sys.stderr.write('usage: %s pdf.lst merge.pdf\n' % __file__)
	quit()

pdf_lst = sys.argv[1]	
pdf_merge = sys.argv[2]	
cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s ' % pdf_merge
for line in open(pdf_lst):
	cmd += '%s ' % line.rstrip('\n')
print(cmd);os.system(cmd)



