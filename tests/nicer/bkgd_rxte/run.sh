#!/bin/sh -f

nipipeline.py \
	--setup_yamlfile $HOPPY_PATH/tests/nicer/bkgd_rxte/setup.yaml \
	--obsid_lstfile $HOPPY_PATH/tests/nicer/bkgd_rxte/obsid.lst
