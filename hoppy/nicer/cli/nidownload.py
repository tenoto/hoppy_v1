#!/usr/bin/env python

__author__ = 'Teru Enoto'
__version__ = '1.03'
# v1.03 : modified for "click" environment 
# v1.02 : modified to use "download_nicer_signleObsID_data.py"
# v1.01 : original version

import sys 
import os 
import click
import hoppy.nicer.nidownload

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="download the NICER target segment summary table.")
@click.option('--username',default=os.getenv('NICERPAGE_USERNAME'),
	help="username to access the NICER team webpage.")
@click.option('--password',default=os.getenv('NICERPAGE_PASSWORD'),
	help="password to access the NICER team webpage.")
@click.option('--outcsvfile',help="output csvfile name.")
def get_target_segment_summary(username,password,outcsvfile):
	nidownloader = hoppy.nicer.nidownload.NicerDownloader()
	nidownloader.get_target_segment_summary(username,password,outcsvfile)

@cli.command(help="download a NICER single ObsID dataset (public-->team).")
@click.argument('obsid')
@click.argument('yyyy_mm')
@click.option('--decrypt',default=os.getenv('NICERDATA_DECRYPT_PASSPHRASE'),
	help="password to decrypt NICER data.")
def get_single_obsid_data(obsid,yyyy_mm,decrypt):
	nidownloader = hoppy.nicer.nidownload.NicerDownloader()
	nidownloader.get_single_obsid_data(obsid,yyyy_mm,decrypt)

@cli.command(help="show a NICER target segment summary table in text format.")
@click.argument('csvfile')
@click.option('--target_name',default=None,help="target name to be shown.")
def show_target_segment_summary(csvfile,target_name):
	nidownloader = hoppy.nicer.nidownload.NicerDownloader()
	nidownloader.show_target_segment_summary(csvfile,target_name=target_name)

@cli.command(help="download multiple NICER ObsID datasets (public-->team).")
@click.option('--username',default=os.getenv('NICERPAGE_USERNAME'),
	help="username to access the NICER team webpage.")
@click.option('--password',default=os.getenv('NICERPAGE_PASSWORD'),
	help="password to access the NICER team webpage.")
@click.option('--decrypt',default=os.getenv('NICERDATA_DECRYPT_PASSPHRASE'),
	help="password to decrypt NICER data.")
def get_multiple_obsid_data(username,password,decrypt):
	nidownloader = hoppy.nicer.nidownload.NicerDownloader()
	nidownloader.get_multiple_obsid_data(username,password,decrypt)

def main():
	cli()

if __name__ == "__main__":
	main()