#!/bin/sh -f 

rm -rf out/multi; mkdir out/multi

$HOPPY_PATH/hoppy/script/make_csv2xspec.py \
	data/filelist2csv.txt \
	data/xspec_multifit_param.yaml \
	out/multi/xspec_multifit.csv

$HOPPY_PATH/hoppy/xspec/xspec_multi_observations.py \
	out/multi/xspec_multifit.csv \
	data/xspec_multifit_param.yaml 




