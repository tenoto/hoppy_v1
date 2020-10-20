#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

df = pd.read_csv("project/scox1/optical/data/vsolj29908_scox1_maehara.dat",
        names=["source","time","magnitude","na"],
        delimiter=" ")
print(df)

magnitude = df.magnitude.str.strip('C').astype(float)

#yyyy = df["time"].astype(str).str[0:4]
#mm = df["time"].astype(str).str[4:6]
#dd = df["time"].astype(str).str[6:8]
#HH = df["time"].astype(str).str[8:10]
#MM = df["time"].astype(str).str[10:12]
#SS = df["time"].astype(str).str[12:14]
#print(yyyy,mm,dd,HH,MM,SS)

#print(df.time)

print(df.time[0:3])
#print(type(df.time[0:3].astype(str)))

#pd.to_datetime(df.time[0:3].astype(str),format="%Y%m%d%H%M%S")

#dtime = pd.to_datetime(df["time"],format="%Y%m%d%H%M%S")
#print(time)

#print(df["time"].astype(str))

print(pd.to_datetime(["20090320"],format="%Y%m%d"))