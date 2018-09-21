#!/bin/sh -f

rm -rf out; rm -f data/J1520-571_00055058g_lc_1day_all_crab.dat

plot_maxi_light_curve.py \
	data/convert_maxilc_to_crabunit.yaml

mv out/J1520-571_00055058g_lc_1day_all_crab.dat data 

plot.py 
