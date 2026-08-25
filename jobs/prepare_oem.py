#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import shutil
import importlib.util

from . import tools
from .tools import oem_regions
from .tools.generate_tracers_xml import generate_tracers_xml
from .tools.ctdas_utilities import (create_prior_all_ones,
                                    create_prior_all_zeros,
                                    create_boundary_prior_all_ones,
                                    create_boundary_prior_separate)

BASIC_PYTHON_JOB = False


def _load_case_script(path):
    """Dynamically load a case-provided Python file as a module.

    Used for custom lambda-/boundary-region generators not covered by
    jobs/tools/oem_regions.py's named strategies: the case's own code gets
    `cfg` passed in directly as a real object, rather than being textually
    substituted the way the remap_icbc bash scripts are (Python source is
    full of literal `{}` -- dict/set literals, f-strings, comprehensions --
    so that .format()-based approach would be far too fragile here).
    """
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _resolve_region_strategy(cfg, config_attr, strategies, default_name,
                             custom_fn_name):
    """Resolve a lambda-/boundary-region strategy from config.yaml.

    `getattr(cfg, config_attr, default_name)` is looked up by name in
    `strategies` first (jobs/tools/oem_regions.py's built-ins); if it isn't
    a known name, it's instead treated as a path (relative to the case
    directory) to a case-provided Python file defining a function called
    `custom_fn_name`.
    """
    name = getattr(cfg, config_attr, default_name)
    if name in strategies:
        return strategies[name]
    module = _load_case_script(cfg.case_path / name)
    return getattr(module, custom_fn_name)


def generate_tracers(cfg):
    """Generate the tracers.xml files (chemtracer definitions for OEM)."""
    tools.create_dir(xml_folder := cfg.case_root / "global_inputs" / "XML",
                     "XML")
    TR_prior = generate_tracers_xml(cfg.tracers,
                                    cfg.CTDAS_nensembles,
                                    restart=False,
                                    propagate_bg=cfg.CTDAS_propagate_bg)
    TR_restart = generate_tracers_xml(cfg.tracers,
                                      cfg.CTDAS_nensembles,
                                      restart=True,
                                      propagate_bg=cfg.CTDAS_propagate_bg)
    with open(xml_folder / "tracers_firstrun.xml", "w",
              encoding="utf-8") as file:
        file.write(TR_prior)
    with open(xml_folder / "tracers_restart.xml", "w",
              encoding="utf-8") as file:
        file.write(TR_restart)
    if cfg.CTDAS_runthrough:
        TR_runthrough_prior = generate_tracers_xml(cfg.tracers,
                                                   cfg.CTDAS_nensembles,
                                                   cfg.CTDAS_nboundaries,
                                                   restart=False,
                                                   runthrough=True)
        TR_runthrough_restart = generate_tracers_xml(cfg.tracers,
                                                     cfg.CTDAS_nensembles,
                                                     cfg.CTDAS_nboundaries,
                                                     restart=True,
                                                     runthrough=True)
        with open(xml_folder / "tracers_runthrough_firstrun.xml",
                  "w",
                  encoding="utf-8") as file:
            file.write(TR_runthrough_prior)
        with open(xml_folder / "tracers_runthrough_restart.xml",
                  "w",
                  encoding="utf-8") as file:
            file.write(TR_runthrough_restart)


def generate_oem_priors(cfg):
    """Generate the OEM lambda regions and initial ensemble priors (all lambdas equal to 1).

    Region generation is a pluggable, named strategy: config.yaml's
    ``prepare_oem.lambda_regions`` / ``prepare_oem.boundary_regions``
    select by name from jobs/tools/oem_regions.py's built-ins --

        per_cell (default), land_ocean_boxes, parent_grid   -- lambda regions
        angular_quadrant (default), compass_walk            -- boundary regions

    -- or, for anything not covered by those, name isn't recognized and is
    instead treated as a path (relative to the case directory) to a
    case-provided Python file defining the same function. See
    jobs/tools/oem_regions.py's docstrings for what each strategy needs
    from cfg.
    """
    tools.create_dir(OEM_folder := cfg.case_root / "global_inputs" / "OEM",
                     "OEM")

    # Interpret lambdas from the YAML file
    lambdas = [
        int(item) for line in cfg.CTDAS_lambdas for item in line.split(',')
    ]

    lambda_regions_fn = _resolve_region_strategy(
        cfg, 'prepare_oem_lambda_regions',
        oem_regions.LAMBDA_REGION_STRATEGIES, 'per_cell',
        'generate_lambda_regions')
    nregs, ncats = lambda_regions_fn(cfg, OEM_folder / "lambdaregions.nc",
                                     lambdas)

    create_prior_all_ones(OEM_folder / "prior_all_ones.nc",
                          nensembles=cfg.CTDAS_nensembles,
                          ncats=max(lambdas),
                          nregs=nregs,
                          propagate_bg=cfg.CTDAS_propagate_bg)
    if cfg.CTDAS_runthrough:
        create_prior_all_zeros(OEM_folder / "prior_all_zeros.nc",
                               nensembles=cfg.CTDAS_nboundaries,
                               ncats=max(lambdas),
                               nregs=nregs)

    # Boundary regions exist to give the ensemble's background (BG)
    # tracers their own spatial perturbation regions -- meaningless (and
    # not worth requiring cdo/iconsub/CTDAS.nboundaries for) unless the
    # case actually has an ensemble tracer to begin with. A plain ICON-ART
    # run with no "-XXX" tracer just skips all of this.
    if any(name.endswith('-XXX') for name in cfg.tracers):
        # Attribution comes from whoever ran the chain, not a fixed name -
        # cfg.user_name/cfg.user_mail are already set from $USER / ~/.forward.
        boundary_regions_fn = _resolve_region_strategy(
            cfg, 'prepare_oem_boundary_regions',
            oem_regions.BOUNDARY_REGION_STRATEGIES, 'angular_quadrant',
            'generate_boundary_regions')
        boundary_regions_fn(cfg, OEM_folder / 'boundary_mask_bg.nc',
                            cfg.CTDAS_nboundaries)

        create_boundary_prior_all_ones(OEM_folder / 'boundary_lambdas_bg.nc',
                                       n_bg_ens=cfg.CTDAS_nboundaries,
                                       nensembles=cfg.CTDAS_nensembles,
                                       propagate_bg=cfg.CTDAS_propagate_bg,
                                       author=cfg.user_name,
                                       email=cfg.user_mail)
        if cfg.CTDAS_runthrough:
            create_boundary_prior_separate(OEM_folder /
                                           'boundary_lambdas_separate.nc',
                                           n_bg_ens=cfg.CTDAS_nboundaries,
                                           author=cfg.user_name,
                                           email=cfg.user_mail)


def main(cfg):
    """
    **Prepare OEM**

    Generates the tracers.xml (chemtracer definitions for OEM) and the OEM
    lambda-region/prior NetCDF inputs. These only need to exist once for the
    whole simulation, so this is a no-op for any chunk after the first.

    Deliberately does *not* set BASIC_PYTHON_JOB, so it runs directly on the
    login node instead of being shipped off to SLURM.

    Parameters
    ----------
    cfg : Config
        Object holding all user-configuration parameters as attributes.
    """
    tools.change_logfile(cfg.logfile)

    if cfg.startdate_sim != cfg.startdate:
        logging.info("Not the first segment, skipping prepare_oem")
        shutil.copy(cfg.logfile, cfg.logfile_finish)
        return

    logging.info("Prepare OEM inputs (tracers.xml, lambda regions, priors)")
    generate_tracers(cfg)
    generate_oem_priors(cfg)
    logging.info("OK")

    # Not a BASIC_PYTHON_JOB, so run_chain.py's own "mark finished" step
    # (shutil.copy(cfg.logfile, cfg.logfile_finish)) never runs for this
    # job -- that only happens for the nested/force_sync invocation
    # BASIC_PYTHON_JOB=True jobs get shipped off to. Do it here instead,
    # so a re-run of this chunk doesn't needlessly redo tracer/OEM setup.
    shutil.copy(cfg.logfile, cfg.logfile_finish)
