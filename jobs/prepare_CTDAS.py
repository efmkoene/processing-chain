#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import xarray as xr
import shutil
import subprocess
from . import tools, prepare_icon
from pathlib import Path  # noqa: F401
from .tools.interpolate_data import create_oh_for_restart, create_oh_for_inicond  # noqa: F401
from .tools.fetch_external_data import fetch_era5, fetch_era5_nudging, fetch_CAMS_CO2, fetch_ICOS_data, fetch_OCO2

BASIC_PYTHON_JOB = False


def main(cfg):
    """
    Prepare CTDAS inversion

    This does the following steps:
    1. Download ERA-5 data for the chosen dates
    2. Download CAMS data for the chosen dates
    3. Interpolate CAMS data to the ERA-5 location (horizontally and vertically)
    4. Download ICOS station data for the chosen dates
    5. Download OCO-2 data for the chosen dates
    6. Prepare the folder output structure

    Parameters
    ----------
    cfg : Config
        Object holding all user-configuration parameters as attributes.
    """
    prepare_icon.set_cfg_variables(cfg)
    print(cfg.print_config())
    tools.change_logfile(cfg.logfile)
    logging.info("Prepare ICON-ART for CTDAS")

    # -- 1. Download ERA5 data and create the initial conditions file
    if cfg.meteo_fetch_era5:
        # -- Fetch ERA5 data
        logging.info(
            f"Times considered now: {cfg.startdate_sim}, {cfg.enddate_sim}, {cfg.CTDAS_step}"
        )
        logging.info("Fetching ERA5 initial data")
        fetch_era5(cfg.startdate_sim, cfg.icon_input_icbc, resolution=0.25)

    # -- 2. Download CAMS CO2 data (for a whole year)
    if cfg.chem_fetch_CAMS:
        fetch_CAMS_CO2(
            cfg.startdate_sim, cfg.icon_input_icbc
        )  # This should be turned into a more central location I think.

    # -- 3. Process data
    # --- ERA5 inicond
    logging.info("Preparing ERA5 preprocessing script for ICON")
    era5_ini_template = cfg.case_path / cfg.meteo_era5_inijob
    era5_ini_job = cfg.icon_input_icbc / cfg.meteo_era5_inijob
    datestr = cfg.startdate_sim.strftime('%Y%m%d%H')
    inicond_filename = cfg.icon_input_icbc / f"era_{datestr}_ini.nc"
    with open(era5_ini_template, 'r') as infile, open(era5_ini_job,
                                                      'w') as outfile:
        outfile.write(infile.read().format(cfg=cfg,
                                           inicond_filename=inicond_filename,
                                           datestr=datestr))
    shutil.copy(cfg.case_path / 'mypartab', cfg.icon_input_icbc / 'mypartab')
    logging.info("Running ERA5 preprocessing script")
    subprocess.run(["bash", era5_ini_job], check=True, stdout=subprocess.PIPE)
    # --- CAMS inicond
    logging.info("Preparing CAMS preprocessing script for ICON")
    cams_ini_template = cfg.case_path / cfg.chem_cams_inijob
    cams_ini_job = cfg.icon_input_icbc / cfg.chem_cams_inijob
    with open(cams_ini_template, 'r') as infile, open(cams_ini_job,
                                                      'w') as outfile:
        outfile.write(infile.read().format(cfg=cfg,
                                           inicond_filename=inicond_filename))
    logging.info("Running CAMS preprocessing script")
    subprocess.run(["bash", cams_ini_job], check=True, stdout=subprocess.PIPE)

    # -- 3. If global nudging, download and process ERA5 and CAMS data
    if cfg.meteo_interpolate_CAMS_to_ERA5:
        for time in tools.iter_hours(cfg.startdate_sim,
                                     cfg.enddate_sim,
                                     step=cfg.meteo_nudging_step):

            # -- Give a name to the nudging file
            timestr = time.strftime('%Y%m%d%H')
            filename = 'era_{timestr}_nudging.nc'.format(timestr=timestr)

            # -- If initial time, copy the initial conditions to be used as boundary conditions
            if time == cfg.startdate_sim:
                shutil.copy(cfg.icon_input_icbc / f'era_{timestr}_ini.nc',
                            cfg.icon_input_icbc / filename)
                continue

            # -- Fetch ERA5 data
            fetch_era5_nudging(time, cfg.icon_input_icbc, resolution=0.25)

            # -- Copy ERA5 processing script (icon_era5_nudging.job) in workdir
            nudging_template = cfg.case_path / cfg.meteo_era5_nudgingjob
            nudging_job = cfg.icon_input_icbc / f'icon_era5_nudging_{timestr}.sh'
            with open(nudging_template, 'r') as infile, open(nudging_job,
                                                             'w') as outfile:
                outfile.write(infile.read().format(cfg=cfg, filename=filename))

            # -- Copy mypartab in workdir
            if not os.path.exists(cfg.case_path / 'mypartab'):
                shutil.copy(cfg.case_path / 'mypartab',
                            cfg.icon_input_icbc / 'mypartab')

            # -- Run ERA5 processing script
            subprocess.run(["bash", nudging_job], check=True, stdout=subprocess.PIPE)

            # -- Copy CAMS processing script (icon_cams_nudging.job) in workdir
            cams_nudging_template = cfg.case_path / cfg.icon_species_nudgingjob
            cams_nudging_job = cfg.icon_input_icbc / f'icon_cams_nudging_{timestr}.sh'
            with open(cams_nudging_template,
                      'r') as infile, open(cams_nudging_job, 'w') as outfile:
                outfile.write(infile.read().format(cfg=cfg, filename=filename))

            # -- Run CAMS processing script
            subprocess.run(["bash", cams_nudging_job], check=True, stdout=subprocess.PIPE)

    # -- 4. Download ICOS CO2 data
    if cfg.obs_fetch_icos:
        # -- This requires you to have accepted the ICOS license in your profile.
        #    So, login to https://cpauth.icos-cp.eu/home/ , check the box, and
        #    copy the cookie token on the bottom as your ICOS_cookie_token.
        fetch_ICOS_data(cookie_token=cfg.ICOS_cookie_token,
                        start_date=cfg.startdate_sim,
                        end_date=cfg.enddate_sim,
                        save_path=cfg.ICOS_path,
                        species=[
                            'co2',
                        ])

    if cfg.obs_fetch_oco2:
        # A user must do the following steps to obtain access to OCO2 data
        # from getpass import getpass
        # import os
        # from subprocess import Popen
        # urs = 'urs.earthdata.nasa.gov'    # Earthdata URL to call for authentication
        # prompts = ['Enter NASA Earthdata Login Username \n(or create an account at urs.earthdata.nasa.gov): ',
        #         'Enter NASA Earthdata Login Password: ']
        # homeDir = os.path.expanduser("~") + os.sep
        # with open(homeDir + '.netrc', 'w') as file:
        #     file.write('machine {} login {} password {}'.format(urs, getpass(prompt=prompts[0]), getpass(prompt=prompts[1])))
        #     file.close()
        # with open(homeDir + '.urs_cookies', 'w') as file:
        #     file.write('')
        #     file.close()
        # with open(homeDir + '.dodsrc', 'w') as file:
        #     file.write('HTTP.COOKIEJAR={}.urs_cookies\n'.format(homeDir))
        #     file.write('HTTP.NETRC={}.netrc'.format(homeDir))
        #     file.close()
        # Popen('chmod og-rw ~/.netrc', shell=True)
        fetch_OCO2(x,
                   y,
                   -8,
                   30,
                   35,
                   65,
                   "/capstor/scratch/cscs/ekoene/temp",
                   product="OCO2_L2_Lite_FP_11.1r")
    logging.info("OK")
