#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

crab_rate = 3.9 

df_maxij1820 = pd.read_csv("project/maxi/lc/data/J1820+071_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_maxij1535 = pd.read_csv("project/maxi/lc/data/J1535-572_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_swiftj0243 = pd.read_csv("project/maxi/lc/data/J0243+614_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_maxij1348 = pd.read_csv("project/maxi/lc/data/J1348-632_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_maxij1631 = pd.read_csv("project/maxi/lc/data/J1631-478_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_scox1 = pd.read_csv("project/maxi/lc/data/J1619-156_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_grs1915 = pd.read_csv("project/maxi/lc/data/J1915+109_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_cygx1 = pd.read_csv("project/maxi/lc/data/J1958+352_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_cygx3 = pd.read_csv("project/maxi/lc/data/J2032+409_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_cirx1 = pd.read_csv("project/maxi/lc/data/J1520-571_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_aqlx1 = pd.read_csv("project/maxi/lc/data/J1911+005_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_4u1608 = pd.read_csv("project/maxi/lc/data/J1612-524_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_4u1634 = pd.read_csv("project/maxi/lc/data/J1634-473_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")
df_4u1705 = pd.read_csv("project/maxi/lc/data/J1708-441_g_lc_1day_all.dat",
	names=["MJD","rate1","error1","rate2","error2","rate3","error3","rate4","error4"],
	delimiter=" ")


plt.style.use("project/maxi/lc/data/supermongo.mplstyle")

mpl.rcParams['axes.linewidth'] = 1.2

fig, axs = plt.subplots(nrows=4,ncols=1,sharex=True,figsize=(10,8))
fig.subplots_adjust(hspace=0)
axs[0].semilogy(df_maxij1535["MJD"], df_maxij1535["rate1"]/crab_rate,'o-',mec='k',label="MAXI J1535-571",markersize=4,markeredgewidth=0)
axs[0].semilogy(df_swiftj0243["MJD"], df_swiftj0243["rate1"]/crab_rate,'o-',mec='k',label="Swift J0243.6+6124",markersize=4,markeredgewidth=0)
axs[0].semilogy(df_maxij1820["MJD"], df_maxij1820["rate1"]/crab_rate,'o-',mec='k',label="MAXI J1820+071",markersize=4,markeredgewidth=0)
axs[0].semilogy(df_maxij1631["MJD"], df_maxij1631["rate1"]/crab_rate,'o-',mec='k',label="MAXI J1631-479",markersize=4,markeredgewidth=0)
axs[0].semilogy(df_maxij1348["MJD"], df_maxij1348["rate1"]/crab_rate,'o-',mec='k',label="MAXI J1348-630",markersize=4,markeredgewidth=0)
axs[0].set_xlim(57754,59142) # 2017-01-01 to 2020-10-20
axs[0].set_ylim(0.03,20)
axs[0].axhline(y=1,linewidth=1,color='k',ls="--")
axs[0].legend(loc=1, title="",framealpha=1.0)
axs[0].set_ylabel(r"X-ray intensity (Crab)") 

axs[1].semilogy(df_scox1["MJD"], df_scox1["rate1"]/crab_rate,'o-',mec='k',label="Sco X-1",markersize=4,markeredgewidth=0)
axs[1].semilogy(df_cygx1["MJD"], df_cygx1["rate1"]/crab_rate,'o-',mec='k',label="Cyg X-1",markersize=4,markeredgewidth=0)
axs[1].semilogy(df_grs1915["MJD"], df_grs1915["rate1"]/crab_rate,'o-',mec='k',label="GRS 1915+105",markersize=4,markeredgewidth=0)
#axs[1].set_ylim(0.1,50)
axs[1].axhline(y=1,linewidth=1,color='k',ls="--")
axs[1].legend(loc=3, title="",framealpha=1.0)
#axs[1].set_ylim(0.01,3)
axs[1].set_ylim(0.005,80)
axs[1].set_ylabel(r"(Crab)") 

axs[2].semilogy(df_cygx3["MJD"], df_cygx3["rate1"]/crab_rate,'o-',mec='k',label="Cyg X-3",markersize=4,markeredgewidth=0)
axs[2].semilogy(df_cirx1["MJD"], df_cirx1["rate1"]/crab_rate,'o-',mec='k',label="Cir X-1",markersize=4,markeredgewidth=0)
axs[2].axhline(y=1,linewidth=1,color='k',ls="--")
axs[2].legend(loc=2, title="",framealpha=1.0)
axs[2].set_ylim(0.005,20)
axs[2].set_ylabel(r"(Crab)") 

axs[3].semilogy(df_aqlx1["MJD"], df_aqlx1["rate1"]/crab_rate,'o-',mec='k',label="Aql X-1",markersize=4,markeredgewidth=0)
axs[3].semilogy(df_4u1608["MJD"], df_4u1608["rate1"]/crab_rate,'o-',mec='k',label="4U 1608-52",markersize=4,markeredgewidth=0)
axs[3].semilogy(df_4u1634["MJD"], df_4u1634["rate1"]/crab_rate,'o-',mec='k',label="4U 1630-472",markersize=4,markeredgewidth=0)
axs[3].set_ylim(0.01,2)
axs[3].axhline(y=1,linewidth=1,color='k',ls="--")
axs[3].legend(loc=2, title="",framealpha=1.0)
axs[3].set_xlabel(r"MJD")
axs[3].set_ylabel(r"(Crab)") 

plt.savefig('maxi_1day_lightcurves.pdf',bbox_inches='tight',transparent=True)








