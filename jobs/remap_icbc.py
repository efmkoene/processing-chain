#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import shutil
import subprocess
from pathlib import Path

from . import tools, prepare_icon

BASIC_PYTHON_JOB = False


def run_bash_script(template, job, **kwargs):
    with job.open('w') as outfile:
        outfile.write(template.read_text().format(**kwargs))
    subprocess.run(["bash", job], check=True)


def main(cfg):
    """
    **Remap ICBC Data**

    Runs the case's ICBC remapping scripts, in order, on the login node.
    Which scripts run -- and what they do (remap ERA5 onto the ICON grid,
    fold in one or more CAMS species, or anything else a case needs) -- is
    entirely up to the case: list them under ``remap_icbc.scripts`` in
    config.yaml, as paths relative to the case directory, e.g.::

        remap_icbc:
            scripts:
                - ICBC/remap_era5.sh
                - ICBC/add_cams_co2.sh

    Each script is rendered with ``.format(cfg=cfg, era5_dir=..., cams_dir=
    ..., outdir=..., workdir=...)`` before being run, so it can use whichever
    of those placeholders (plus any ``cfg.xxx`` attribute) it needs -- e.g. a
    case fetching several CAMS species could list one script per species, or
    a single script handling all of them.

    Parameters
    ----------
    cfg : Config
        Object holding all user-configuration parameters as attributes.
    """
    prepare_icon.set_cfg_variables(cfg)
    tools.change_logfile(cfg.logfile)
    logging.info("Remap ICBC data")

    scripts = getattr(cfg, 'remap_icbc_scripts', [])
    if not scripts:
        logging.info("No remap_icbc.scripts configured, nothing to do")
        shutil.copy(cfg.logfile, cfg.logfile_finish)
        return

    workdir = cfg.chain_root / "remap_icbc_work"
    tools.create_dir(workdir, "remap_icbc work directory")

    kwargs = dict(
        cfg=cfg,
        era5_dir=cfg.case_root / "global_inputs" / "ERA5",
        cams_dir=cfg.case_root / "global_inputs" / "CAMS",
        outdir=cfg.icon_input_icbc,
        workdir=workdir,
    )

    for script_rel in scripts:
        template = cfg.case_path / script_rel
        job = workdir / Path(script_rel).name
        logging.info(f"Running {template}")
        run_bash_script(template, job, **kwargs)

    logging.info("OK")

    # Not a BASIC_PYTHON_JOB, so run_chain.py's own "mark finished" step
    # (shutil.copy(cfg.logfile, cfg.logfile_finish)) never runs for this
    # job -- that only happens for the nested/force_sync invocation
    # BASIC_PYTHON_JOB=True jobs get shipped off to. Do it here instead,
    # so a re-run of this chunk doesn't needlessly redo the remap.
    shutil.copy(cfg.logfile, cfg.logfile_finish)
