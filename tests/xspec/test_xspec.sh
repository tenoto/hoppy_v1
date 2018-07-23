#!/bin/sh -f 

rm -rf out;
python ~/work/drbv1/soft/git/hoppy/hoppy/xspec/xspec.py \
	data/ni1200250101_o1u20cut_clscr_gtisel.pha \
	-o out \
	-b data/ni1200250101_o1u20cut_BGMod_3C50.pha \
	-m script/input/phabs_pegpwrlw.xcm \
	-s 5 -n 80 --fitemin 0.3 --fitemax 10.0 \
	--rateband 0.8-6.0,2.0-10.0 \
	--fluxband 0.8-6.0,2.0-10.0	
