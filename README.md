# HOPPY: High-energy Observatory Pipelines via PYthon

###### tags: `python` `library`

[![hackmd-github-sync-badge](https://hackmd.io/KktKqjMbT3Gj3zNO0Jpk6w/badge)](https://hackmd.io/KktKqjMbT3Gj3zNO0Jpk6w)

:beer::beer::beer: https://github.com/tenoto/hoppy
 
このパッケージは、Ｘ線データ解析でよく使う処理を自動化するためのライブラリやスクリプトの集合体です。

## 初期設定

### Python 環境

`pyenv` と `pyenv-virtualenv` の環境で動かせるのが理想です。

### パス設定

シェル設定ファイル (e.g., ~/.zshrc) に以下の

```
export HOPPY_PATH="Path/To/User/Hoppy/Directory"
alias hoppyinit="source $HOPPY_PATH/setenv/setenv.bashrc"
```

毎回、コマンドプロンプトを立ち上げるたびに、

```
heainit
hoppyinit
```

とすれば動かせるようになるはずです。Python の必要なライブラリは要求されたら、`pip install xxx` で随時入れてください。

(本当は python ライブラリのバージョン管理して、環境ごとダウンロードできるようにしたいのですが、時間が取れないので、とりあえず...）

## Xspec 

Xspec の自動解析を、コマンドラインをそのまま書き下して実装。

### ライブラリ

- `xspec.py`:

### コマンド (CLI)

- `get_xspec_rate.py`: Xspec の `show rate` コマンドで表示されるレートを読んでくれる。
- `xspec_fit.py`:　ひとつの Xspec のフィッティングを実行する。
- `make_csv2xspec.py`: 入力するスペクトルファイル、rmf, arf のレスポンスファイル、パラメータの yaml ファイルを入力したら、`xspec_multi_fit.py` の入力ファイルに変換してくれる。
- `xspec_multi_fit.py`: `xspec_fit.py` をシリーズで複数回実行する。

## NICER 

```
hoppy/nicer
├── cli (コマンドライン・スクリプト群)
│   ├── niwget.py
│   ├── ...
│   └── niauto.py
├── nicer.py (モジュール)
```

### ライブラリ

- `nicer.py`: メインのライブラリ。`NicerElf` クラスは、寝ている間に仕事をしてくれる妖精さん(Elf)で、複数の ObsID のデータを扱う。`NicerObsID` と `NicerGTI` はそれぞれ、ObsID 単位、GTI 単位でのデータを扱う。

### コマンド (CLI)

#### ダウンロード関係
- `niget_target_summary_sheet.py`: NICER Team のターゲット天体と ObsID リスト(Target Summary Sheet)をダウンロードします。チームの ID と Password が必要になります。
- `nishow_target_segment_sheet.py`: Target Summary Sheet から、該当する天体もしくは ObsID の情報を抜き出して表示してくれます。
- `niget_yyyymm.py`: ObsID を指定し、該当するデータが保存されているディレクトリ（観測年と月）を教えてくれる。
- `niwget.py`: 指定した ObsID もしくは、ターゲット名のデータをすべてダウンロードしてくる。HEASARC の公式アーカイブを探し、もしなければ NICER チーム内のアーカイブを確認する。

#### パイプライン処理

- `nipipeline.py`: 設定条件ファイルのレシピにしたがって、複数の ObsID についてパイプライン処理をしていきます。簡易的にスペクトル生成の結果まで作りますが、フィットは別にしてあります。
- `niauto.py`: データダウンロード、パイプライン処理、スペクトル・フィットまでを自動化してある。

## 時間変動解析 (未実装・実装中)

### ライブラリ

- `timeseries.py`: 

### コマンド (CLI)

- `tmsearchburst.py`: イベント形式のデータで、一定幅の時間変動を示すライトカーブから、それに当てはまらな増大を検出する。

## Reference
- Python package lesson https://github.com/BillMills/pythonPackageLesson (see demo python package)
- PEP 8 -- Style Guide for Python Code https://www.python.org/dev/peps/pep-0008/?

## History
- 2020-08-17 modified (Teruaki Enoto)
- 2018-11-14 new version 0.1.0 is created for CLI interface using click.



<!--
This package includes libararies and scripts for X-ray analyses of High-energy Observatory Pipeines of PYthon (HOPPY). 

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


## NICER 

```
hoppy/nicer
├── cli
│   ├── nibarytime.py
│   ├── niget_target_segment_sheet.py
│   ├── niget_yyyymm.py
│   ├── nipipe01_nicerl2.py
│   ├── nishow_target_segment_sheet.py
│   ├── nitimeconv.py
│   ├── niwget.py
│   └── plot_mitbgd3c50_photfile.py
├── nicer.py
├── nievent.py
```

* nicer.py: main library of the NICER data analysis framework. The "NicerElf" class is organizign the framwork. There are two main calsses NicerObsID and NicerGTI. 
* nipipeline.py: command line to run the nicer.py process.
* niwget.py
* niget_yyyymm.py
* nishow_target_segment_sheet.py
* niget_target_segment_sheet.py 

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
 
-->
