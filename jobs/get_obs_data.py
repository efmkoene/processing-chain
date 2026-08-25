#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import shutil
from datetime import timedelta

from . import tools
from .tools.fetch_obs_data import (fetch_ICOS_data, process_ICOS_data,
                                   fetch_OCO2_data, process_OCO2_data,
                                   load_icon_polygon_erode_buffer)

BASIC_PYTHON_JOB = False


def main(cfg):
    """
    **Fetch Obs Data**

    Downloads ICOS station and OCO-2 satellite CO2 observations for the
    current chunk, masked/matched against the ICON domain grid.

    Deliberately does *not* set BASIC_PYTHON_JOB, so it runs directly on the
    login node instead of being shipped off to SLURM: both fetches spend
    most of their time waiting on a remote service, not doing local compute.

    Parameters
    ----------
    cfg : Config
        Object holding all user-configuration parameters as attributes.
    """
    tools.change_logfile(cfg.logfile)
    logging.info("Fetch obs data (ICOS / OCO-2)")

    fetch_start = cfg.startdate_sim
    fetch_end = cfg.enddate_sim + timedelta(days=1)
    grid_file = cfg.input_files_dynamics_grid_filename
    n_layers = getattr(cfg, 'CTDAS_obs_n_layers', 12)

    if getattr(cfg, 'CTDAS_obs_ICOS_fetch', False):
        fetch_ICOS_data(start_date=fetch_start.strftime("%d-%m-%Y"),
                        end_date=fetch_end.strftime("%d-%m-%Y"),
                        save_path=cfg.CTDAS_obs_ICOS_path,
                        species=['co2'])

        icos_dir = cfg.case_root / "global_inputs" / "ICOS"
        tools.create_dir(icos_dir, "ICOS input files")
        icon_poly = load_icon_polygon_erode_buffer(grid_file, n_layers=n_layers)
        process_ICOS_data(ICOS_obs_folder=cfg.CTDAS_obs_ICOS_path,
                          start_date=fetch_start,
                          end_date=fetch_end,
                          output_folder=icos_dir,
                          icon_poly=icon_poly)

    if getattr(cfg, 'CTDAS_obs_OCO2_fetch', False):
        fetch_OCO2_data(DATE_WINDOW=(fetch_start.strftime("%Y-%m-%d"),
                                    fetch_end.strftime("%Y-%m-%d")),
                        OUT_DIR=cfg.CTDAS_obs_OCO2_path,
                        ICON_GRID_PATH=grid_file,
                        n_layers=n_layers)

        oco2_dir = cfg.case_root / "global_inputs" / "OCO2"
        tools.create_dir(oco2_dir, "OCO-2 output")
        process_OCO2_data(OCO2_obs_folder=cfg.CTDAS_obs_OCO2_path,
                          ICON_grid_file=grid_file,
                          start_date=fetch_start,
                          end_date=fetch_end,
                          output_folder=oco2_dir)

    logging.info("OK")

    # Not a BASIC_PYTHON_JOB, so run_chain.py's own "mark finished" step
    # (shutil.copy(cfg.logfile, cfg.logfile_finish)) never runs for this
    # job -- that only happens for the nested/force_sync invocation
    # BASIC_PYTHON_JOB=True jobs get shipped off to. Do it here instead,
    # so a re-run of this chunk doesn't needlessly re-fetch.
    shutil.copy(cfg.logfile, cfg.logfile_finish)
