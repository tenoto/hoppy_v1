#!/bin/sh -f 

rm -rf out_multi; mkdir out_multi

python ../../../hoppy/script/make_csv2xspec.py \
	data/filelist2csv.txt \
	data/xspec_multifit_param.yaml \
	out_multi/xspec_multifit.csv

#python ../../../hoppy/xspec/xspec_multi_observations.py \
#	data/filelist2csv.txt \
#	data/xspec_multifit_param.yaml 




