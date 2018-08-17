#!/bin/sh -f 

# heasarc public data
download_nicer_singleObsID_data.py 1012010110 2017_09

# team internal data 
download_nicer_singleObsID_data.py 1070010166 2018_02 --gpg_decryption_password "XXXX"

# error 
download_nicer_singleObsID_data.py 1012010110 2017_08
