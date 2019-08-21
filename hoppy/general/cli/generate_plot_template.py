#!/usr/bin/env python

import os 
import sys

if len(sys.argv) != 2:
	print("error %s file.txt" % sys.argv[0])
	quit()

infile = sys.argv[1]
basename = os.path.splitext(os.path.basename(infile))[0]
script_file = "plot_%s.py" % basename 
outpdf = "%s.pdf" % basename

f = open(script_file,'w')
dump = """#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
"""
dump += "hdu = pd.read_csv('%s',header=0,delim_whitespace=True)" % infile
dump += """
plt.style.use('https://raw.githubusercontent.com/tenoto/repository/master/style/matplotlib/aspect_8to6.mplstyle')

fig, ax = plt.subplots()
fig.patch.set_alpha(0.0)
plt.tick_params(labelsize=18,direction='in')
plt.errorbar(hdu["X"],hdu["Y"],
	xerr=hdu["Xerr"],yerr=hdu["Yerr"],
	marker='.',ls='',drawstyle='steps-mid',color='C3')
ax.set_xlabel(r"Energy (keV)")
ax.set_ylabel(r"Counts s$^{-1}$ keV$^{-1}$ cm$^{-2}$") 
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.2,15.0)
ax.set_ylim(2e-5,0.02)
ax.patch.set_alpha(0.0)
"""
dump += "plt.savefig('%s',bbox_inches='tight')" % outpdf
f.write(dump)
f.close()

os.system('chmod +x %s' % script_file)

