#!/usr/bin/env python

import os 
import sys 
if len(sys.argv) != 2:
	sys.stderr.write('usage: %s file.ps\n' % __file__)
	quit()

psfile = sys.argv[1]	
if os.path.splitext(psfile)[-1] != '.ps':
	sys.stderr.write('input file is not ps file\n' % psfile)
	quit()
cmd = 'ps2pdf %s' % psfile
print(cmd);os.system(cmd)

pdffile = os.path.basename(psfile).replace('.ps','.pdf')
cmd = 'mv %s %s; rm -f %s' % (pdffile,os.path.dirname(psfile),psfile)
print(cmd);os.system(cmd)


