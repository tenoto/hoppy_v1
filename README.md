HOPPY (High-energy Observatory Pipelines via PYthon)
===
This package includes libararies and scripts for X-ray analyses of High-energy Observatory Pipelines of PYthon (HOPPY). 

https://github.com/tenoto/hoppy

Latest version 0.1.0 

# Setup 
Write following lines to the initialization setup (e.g., ~/.zshrc).

```
export HOPPY_PATH="/Users/enoto/work/drbv1/soft/git/hoppy"
alias hoppyinit="cd $HOPPY_PATH; source setenv/setenv.bashrc; pipenv shell"
```

Everytime, you need to run the following command line inputs, 

```
heainit
hoppyinit
```


## Files and Libraries
* hoppy
    * nicer 
        * [nievent.py](https://github.com/tenoto/hoppy/blob/master/hoppy/nicer/nievent.py) library 
    * maxi 
    * nustar
    * physics
    * plot
    * rxte 
    * script
    * swift
    * timing 
    * xmm
    * xspec
    * xte
* gallery (example of plots and illustrations)
    * QDP
    * 2018
* setenv/setenv.bashrc (environmental setups)
* tests (test scripts)

## Structure

```
hoppy
├── LICENSE
├── MANIFEST.in
├── README.md
├── dist
├── hoppy
│   ├── general
│   │   └── __init__.py
│   ├── nicer
│   │   └── __init__.py
│   ├── plot
│   │   └── __init__.py
│   └── swift
│       └── __init__.py
└── setup.py
```

## Class Description

- XrayObservationDatabase
    - attribute
        - name
        - directory_path
        - csvfile
        - htmlfile
    - methods
        - run_pipeline
- XrayObservation(XrayObservationDatabase)
    - attribute
        - target 
        - satellite
        - obsid 
        - name (target_satellite_obsid)
        - directory_path (target/satellite/obsid)
        - yamlfile
        - detector
        - obsmode
        - start_mjd
        - stop_mjd
        - start_yyyymmdd
        - start_hhmmss
        - stop_yyyymmdd
        - stop_hhmmss
        - exposure
        - ra_obs
        - dec_obs
        - rate
        - rate_error
        - source_count
        - backgrnd_count
        - soruce_pha
        - backgrnd_pha
        - rmffile
        - arffile
        - grppha_min_significance
        - grppha_max_bins
        - grppha_emin
        - grppha_emax
        - flag_generate_spectra (done,error,none)
        - flag_fit_spectra (done,error,none)
        - xspec_model_xcmfile
        - parnum_fixed
        - parnum_error
        - parameters ([parnum,value,error_min,error_max]...])
    - methods
        - generate_spectra
        - fit_spectrum

## Reference
- Python package lesson https://github.com/BillMills/pythonPackageLesson (see demo python package)
- PEP 8 -- Style Guide for Python Code https://www.python.org/dev/peps/pep-0008/?

## History
- 2018-11-14 new version 0.1.0 is created for CLI interface using click.
