#!/bin/sh -f

nibgdlc3c50.py \
	--bkgidxdir '/Users/enoto/work/soft/heasoft/nicer/nibackgen3c50/develop/202007_RGv6/bg_models_3C50' \
	--bkglibdir '/Users/enoto/work/soft/heasoft/nicer/nibackgen3c50/develop/202007_RGv6/bg_models_3C50' \
	--gainepoch 2020 \
	--ufafile 'out/bkgd_rxte/1012020134/xti/event_cl/ni1012020134_0mpu7_ufa.evt' \
	--clfile 'out/bkgd_rxte/1012020134/xti/event_cl/ni1012020134_0mpu7_cl.evt' \
	--dtmin 20.0 \
	--dtmax 60.0 \
	--hbgcut 0.5 \
	--fpmofflist 14,34 \
	--tbin 60
