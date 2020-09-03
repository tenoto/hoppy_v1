#!/bin/sh -f

# Crab Pulsar
# 2018-04-07
# exposure 2945.76

 # number of events: 1e+5 counts
 # start time : MET = 134559840.000061
 # stop time : MET = 134559849.016345
 # effective exposure = ~9 sec

xselect<<EOF
xsel
read event  1013010131/xti/event_cl/ni1013010131_0mpu7_cl.evt.gz ./
yes
filer time scc 
134559840.0,134559860
x
extract event 
save event
ni1013010131_0mpu7_cl_20sec.evt
yes
extract spec
save spec
ni1013010131_0mpu7_cl_20sec.pha
exit
no
EOF

barycorr infile=ni1013010131_0mpu7_cl_20sec.evt \
outfile=ni1013010131_0mpu7_cl_20sec_bary.evt \
orbitfiles="1013010131/auxil/ni1013010131.orb.gz" \
ra=83.633218 dec=22.014464 refframe=ICRS ephem=JPLEPH.430

powspec
ni1013010131_0mpu7_cl_20sec_bary.evt
-
1e-6
2000000
10
0
ni1013010131_0mpu7_cl_20sec_bary.pow
yes
/xw
