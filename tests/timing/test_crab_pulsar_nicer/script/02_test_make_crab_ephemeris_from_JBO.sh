#!/bin/sh -f

make_crab_ephemeris.py 57980

mv crab_JBO_ephemeris_MJD57980.par data

rm -f crab2.txt 