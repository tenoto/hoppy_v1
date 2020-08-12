#!/bin/sh -f 

rm -rf out/multi; mkdir -p out/multi

$HOPPY_PATH/hoppy/xspec/cli/make_csv2xspec.py \
	data/filelist2csv.txt \
	data/xspec_multifit_param.yaml \
	out/multi/xspec_multifit.csv

$HOPPY_PATH/hoppy/xspec/cli/xspec_multi_fit.py \
	out/multi/xspec_multifit.csv \
	data/xspec_multifit_param.yaml 




