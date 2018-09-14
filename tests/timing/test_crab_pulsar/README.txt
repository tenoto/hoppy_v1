
Crab observation [X-ray and radio]

Date : 2017-08-09 (MJD 57974, doy 221)
NICER ObsID : 1013010104 (Exposure 10.1 ks, UTC 17:13-23:53)
Usuda radio observatory: UTC 16:45-23:59 UTC

-------------------------
Extract sample event file
-------------------------
Usuda 2nd GTI from 2017221_UsdSch0123_GTItbl.txt
				Start 				End
Usuda-MJD-UTC? 	57974.750000000000  57974.791655080973 (duraiton 3599.286457)
NICER-MET 		113767202.000		113770800.999

fselect ni1013010104_0mpu7_cl.evt.gz \
	ni1013010104_0mpu7_cl_crab_sample.evt \
	"(TIME >= 113767202.000) && (TIME < 113770800.999)"

Number of the original (ni1013010104_0mpu7_cl.evt.gz) : 93,085,846 (1.9 G)
Number of the filtered (ni1013010104_0mpu7_cl_crab_sample.evt) : 12,118,944 (451 M)

further filtering for 10-minute duration 

fselect ni1013010104_0mpu7_cl.evt.gz \
	ni1013010104_0mpu7_cl_crab_10min.evt \
	"(TIME >= 113769700.0) && (TIME < 113770300.0)"

Number of the filtered (ni1013010104_0mpu7_cl_crab_10min.evt) : 6,517,599 (243 M)

fselect ni1013010104_0mpu7_cl.evt.gz \
	ni1013010104_0mpu7_cl_crab_1min.evt \
	"(TIME >= 113769700.0) && (TIME < 113769760.0)"

Number of the filtered (ni1013010104_0mpu7_cl_crab_1min.evt) : 653,128 (24 M)

---------------------------------
(A.1) JBO monthly ephemeris (DE200)
---------------------------------
./test_make_crab_ephemeris_from_JBO.sh

photonphase \
	--barytime --absphase --ephem DE200 \
	--orbfile data/ni1013010104.orb \
	--outfile ni1013010104_0mpu7_cl_crab_1min_absphase_JBO.evt \
	--plotfile ni1013010104_0mpu7_cl_crab_1min_absphase_JBO.png \
	ni1013010104_0mpu7_cl_crab_1min.evt \
	crab_JBO_ephemeris_MJD57980.par 
450.94s user 9.63s system 98% cpu 7:45.35 total

---------------------------------
(A.2) Barycorr with the same parameter as A.1
---------------------------------
barycorr \
	infile=ni1013010104_0mpu7_cl_crab_1min.evt \
	outfile=ni1013010104_0mpu7_cl_crab_1min_barycorr_JBO.evt \
	ra=83.633218 dec=22.014464 \
	orbitfiles=data/ni1013010104.orb \
	refframe=FK5 ephem=JPLEPH.200 

First event : 
TIME (before barycorr) : 113769700.000001 
photonphase : BARY_TIME : 57974.7762127998
			: ABS_PHASE : -13377348
			: PHOTON_PHASE : 0.481783664530667
barycorr : TIME (after barycorr) : 113769397.601902			
		==> 57974.77541206 (converted from xTime)
		==> 57974.7762127998

MJDREFI + MJDREFF + TIME/86400.0. 	

(before barycorr)
MJDREFI = 56658	
MJDREFF = 0.000777592592592593

(after barycorr)
MJDREFI = 56658	
MJDREFF = 0.000777592592592593 

---------------------------------
(B) Usuda (DE430)
---------------------------------
barycorr infile=data/ni1013010104_0mpu7_cl_crab_1min.evt \
	outfile=ni1013010104_0mpu7_cl_crab_1min_bary.evt \
	ra=83.633218 dec=22.014464 \
	orbitfiles=data/ni1013010104.orb \
	refframe=ICRS ephem=JPLEPH.430 

faddphase_nu.py \
	-i ni1013010104_0mpu7_cl_crab_1min_bary.evt \
	-o ni1013010104_0mpu7_cl_crab_1min_bary_phase.evt \
	-n 29.639601201518 \
	-d -3.687105e-10 \
	-e 113702332.821496

fplot_pulseprofile.py \
 	ni1013010104_0mpu7_cl_crab_1min_bary_phase.evt \
 	--outfits ni1013010104_0mpu7_cl_crab_1min_bary_phase_pls.fits \
 	--nbin 100 --colname PHASE	



















