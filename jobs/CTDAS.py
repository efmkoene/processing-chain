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
from .tools.fetch_external_data import fetch_era5, fetch_era5_nudging

BASIC_PYTHON_JOB = False


def main(cfg):
    """
    Prepare CTDAS inversion

    This does the following steps:
    1. Run the first day (spin-up)
    2. Start CTDAS

    Parameters
    ----------
    cfg : Config
        Object holding all user-configuration parameters as attributes.
    """
    prepare_icon.set_cfg_variables(cfg)
    tools.change_logfile(cfg.logfile)
    logging.info("Prepare ICON-ART for global simulations")

    # -- Download ERA5 data and create the inicond file
    if cfg.era5_inicond and cfg.lrestart == '.FALSE.':
        # -- Fetch ERA5 data
        fetch_era5(cfg.startdate_sim, cfg.icon_input_icbc)

        # -- Copy ERA5 processing script (icon_era5_inicond.job) in workdir
        with open(cfg.icon_era5_inijob) as input_file:
            to_write = input_file.read()
        output_file = os.path.join(cfg.icon_input_icbc, 'icon_era5_inicond.sh')
        with open(output_file, "w") as outf:
            outf.write(to_write.format(cfg=cfg))
