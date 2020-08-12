#!/bin/sh -f

python $HOPPY_PATH/hoppy/xspec/cli/get_xspec_rate.py \
	data/ni1200250101_o1u20cut_clscr_gtisel.pha 1.0 8.0 \
	-b data/ni1200250101_o1u20cut_BGMod_3C50.pha \
	-r data/nicer_v1.02.rmf \
	-a data/ni_xrcall_onaxis_v1.02_wo_detid14_34.arf \

