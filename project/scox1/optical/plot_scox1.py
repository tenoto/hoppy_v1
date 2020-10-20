#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates

df = pd.read_csv("project/scox1/optical/data/vsolj29908_scox1_maehara.dat",
        names=["source","time","magnitude","na"],
        delimiter=" ")
print(df)

magnitude = df.magnitude.str.strip('C').astype(float)

yyyy = df["time"].astype(str).str[0:4]
mm = df["time"].astype(str).str[4:6]
dd = df["time"].astype(str).str[6:8]
HH = (df["time"].astype(str).str[8:10].astype(int) - 9).astype(str)
MM = df["time"].astype(str).str[10:12]
SS = df["time"].astype(str).str[12:14]

str_date = yyyy + mm + dd + HH + MM + SS

obstime = pd.to_datetime(str_date,format="%Y%m%d%H%M%S")

fig, ax = plt.subplots(nrows=1,ncols=1,sharex=True,figsize=(12,4))
plt.gca().invert_yaxis()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.step(obstime, magnitude)
ax.set_xlabel(r"UT Time on March 20, 2009") 
ax.set_ylabel(r"Magnitude") 
#ax.plot(obstime, magnitude,'o-',mec='k',label="Sco X-1",markersize=4,markeredgewidth=1)
#axs[0].set_xlim(57754,59142) # 2017-01-01 to 2020-10-20
#axs[0].set_ylim(0.03,20)
#axs[0].axhline(y=1,linewidth=1,color='k',ls="--")


plt.savefig('scox1_optical_200903.pdf',bbox_inches='tight',transparent=True)








