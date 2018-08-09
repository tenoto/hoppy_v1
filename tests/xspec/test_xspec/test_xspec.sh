#!/bin/sh -f 

rm -rf out;
python $HOPPY_PATH/hoppy/xspec/xspec.py \
	data/ni1200250101_o1u20cut_clscr_gtisel.pha \
	-o out/ni1200250101 \
	-b data/ni1200250101_o1u20cut_BGMod_3C50.pha \
	-r data/nicer_v1.02.rmf \
	-a data/ni_xrcall_onaxis_v1.02_wo_detid14_34.arf \
	-m data/phabs_pegpwrlw.xcm \
	-s 5 -n 80 --fitemin 0.3 --fitemax 10.0 \
	--rateband 0.8-6.0,2.0-10.0 \
	--fluxband 0.8-6.0,2.0-10.0	
