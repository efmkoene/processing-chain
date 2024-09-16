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
from .tools.fetch_external_data import fetch_era5, fetch_era5_nudging, fetch_CAMS_CO2, fetch_ICOS

BASIC_PYTHON_JOB = True


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
    tools.change_logfile(cfg.logfile)
    logging.info("Prepare ICON-ART for CTDAS")

    # -- 1. Download ERA5 data and create the initial conditions file
    if cfg.era5_inicond:
        # -- Fetch ERA5 data
        fetch_era5(cfg.startdate_sim, cfg.icon_input_icbc)

        # -- Copy ERA5 processing script (icon_era5_inicond.job) in workdir
        with open(cfg.icon_era5_inijob) as input_file:
            to_write = input_file.read()
        output_file = os.path.join(cfg.icon_input_icbc, 'icon_era5_inicond.sh')
        with open(output_file, "w") as outf:
            outf.write(to_write.format(cfg=cfg))

        # -- Copy mypartab in workdir
        shutil.copy(
            os.path.join(os.path.dirname(cfg.icon_era5_inijob), 'mypartab'),
            os.path.join(cfg.icon_input_icbc, 'mypartab'))

        # -- Run ERA5 processing script
        process = subprocess.Popen([
            "bash",
            os.path.join(cfg.icon_input_icbc, 'icon_era5_inicond.sh')
        ],
                                   stdout=subprocess.PIPE)
        process.communicate()

    # -- 2. Download CAMS CO2 data
    if cfg.cams_inicond:
        fetch_CAMS_CO2(cfg.startdate_sim, cfg.icon_input_icbc)

    # ((( Is there an interpolation step missing for inicond? I assume so...)))

    # -- 3. If global nudging, download and process ERA5 and CAMS data
    if cfg.era5_cams_nudging:
        for time in tools.iter_hours(cfg.startdate_sim,
                                     cfg.enddate_sim,
                                     step=cfg.nudging_step):

            # -- Give a name to the nudging file
            timestr = time.strftime('%Y%m%d%H')
            filename = 'era_{timestr}_nudging.nc'.format(timestr=timestr)

            # -- If initial time, copy the initial conditions to be used as boundary conditions
            if time == cfg.startdate_sim and cfg.era5_inicond:
                shutil.copy(cfg.input_files_scratch_inicond_filename,
                            os.path.join(cfg.icon_input_icbc, filename))
                continue

            # -- Fetch ERA5 data
            fetch_era5_nudging(time, cfg.icon_input_icbc)

            # -- Copy ERA5 processing script (icon_era5_nudging.job) in workdir
            with open(cfg.icon_era5_nudgingjob) as input_file:
                to_write = input_file.read()
            output_file = os.path.join(
                cfg.icon_input_icbc, 'icon_era5_nudging_{}.sh'.format(timestr))
            with open(output_file, "w") as outf:
                outf.write(to_write.format(cfg=cfg, filename=filename))

            # -- Copy mypartab in workdir
            if not os.path.exists(os.path.join(cfg.icon_input_icbc,
                                               'mypartab')):
                shutil.copy(
                    os.path.join(os.path.dirname(cfg.icon_era5_nudgingjob),
                                 'mypartab'),
                    os.path.join(cfg.icon_input_icbc, 'mypartab'))

            # -- Run ERA5 processing script
            process = subprocess.Popen([
                "bash",
                os.path.join(cfg.icon_input_icbc,
                             'icon_era5_nudging_{}.sh'.format(timestr))
            ],
                                       stdout=subprocess.PIPE)
            process.communicate()

            # -- Copy CAMS processing script (icon_cams_nudging.job) in workdir
            with open(cfg.icon_species_nudgingjob) as input_file:
                to_write = input_file.read()
            output_file = os.path.join(
                cfg.icon_input_icbc, 'icon_cams_nudging_{}.sh'.format(timestr))
            with open(output_file, "w") as outf:
                outf.write(to_write.format(cfg=cfg, filename=filename))

            # -- Run CAMS processing script
            process = subprocess.Popen([
                "bash",
                os.path.join(cfg.icon_input_icbc,
                             'icon_cams_nudging_{}.sh'.format(timestr))
            ],
                                       stdout=subprocess.PIPE)
            process.communicate()

    # -- 4. Download ICOS CO2 data
    if cfg.fetch_ICOS:
        # -- This requires you to have accepted the ICOS license in your profile.
        #    So, login to https://cpauth.icos-cp.eu/home/ , check the box, and
        #    copy the cookie token on the bottom as your ICOS_cookie_token.
        fetch_ICOS(cookie_token=cfg.ICOS_cookie_token,
                   start_date=cfg.startdate_sim,
                   end_date=cfg.enddate_sim,
                   save_path=cfg.ICOS_path,
                   species=[
                       'co2',
                   ])

    logging.info("OK")
