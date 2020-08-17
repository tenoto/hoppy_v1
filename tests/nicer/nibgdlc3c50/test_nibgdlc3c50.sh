#!/bin/sh -f

nibgdlc3c50.py \
	--clfile 'out/bkgd_rxte/1012040184/xti/event_cl/ni1012040184_0mpu7_cl.evt' \
	--ufafile 'out/bkgd_rxte/1012040184/xti/event_cl/ni1012040184_0mpu7_ufa.evt' \
	--tbinfile 'ni1012040184_0mpu7_cl_tbin.csv' \
	--bkgidxdir '/Users/enoto/work/soft/heasoft/nicer/nibackgen3c50/develop/202007_RGv6/bg_models_3C50' \
	--bkglibdir '/Users/enoto/work/soft/heasoft/nicer/nibackgen3c50/develop/202007_RGv6/bg_models_3C50' \
	--gainepoch 2020 \
	--dtmin 20.0 \
	--dtmax 60.0 \
	--hbgcut 0.1 \
	--fpmofflist 14,34 \
	--rmffile 'resp/nicer-rmf6s-teamonly-array50_wo14_34.rmf' \
	--arffile 'resp/nicer-consim135p-teamonly-array50_wo14_34.arf' \
	--ratebands 0.3-2.0,0.0-0.3,2.0-8.0,8.0-12.0 \
	--outcsvfile 'ni1012040184_0mpu7_cl_nibgdlc3c50.csv' 

