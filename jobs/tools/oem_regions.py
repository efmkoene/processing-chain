"""Named strategies for OEM lambda-region / boundary-region generation.

prepare_oem.py selects one of these by name via config.yaml
(``prepare_oem.lambda_regions`` / ``prepare_oem.boundary_regions``); a name
not found here is instead treated as a path to a case-provided Python file
defining the same function, for anything not covered by the built-ins below.

Every function takes ``cfg`` directly (never templated into source text --
Python is full of literal ``{}``, so the .format()-based approach used for
the bash ICBC-remap scripts would be far too fragile here) and follows one
of two signatures:

    generate_lambda_regions(cfg, output_path, lambdas) -> (nregs, ncats)
    generate_boundary_regions(cfg, output_path, n_bg_ens) -> None
"""
import subprocess

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree


def lambda_regions_per_cell(cfg, output_path, lambdas):
    """Every grid cell is its own region -- the finest possible resolution,
    and the default. Needs cfg.input_files_dynamics_grid_filename.
    """
    ds = xr.open_dataset(cfg.input_files_dynamics_grid_filename)
    ncells = ds.cell.size
    nregs = ncells
    categories = np.arange(1, len(lambdas) + 1)
    regions = np.arange(nregs)
    cells = np.arange(ncells) + 1

    ds_cells = xr.Dataset(data_vars={
        'REG': (['cell'], regions),
        'Lambda_indicies': (['cat'], lambdas)
    },
                          coords={
                              'cell': (['cell'], cells),
                              'cat': (['cat'], categories)
                          },
                          attrs={'author': 'Processing Chain'})
    try:
        ds_cells.to_netcdf(output_path,
                           encoding={
                               'REG': {
                                   'dtype': 'int32'
                               },
                               'cell': {
                                   'dtype': 'int32'
                               }
                           })
    except:
        print("File currently open. Please close the file and try again.")
    print(f"Lambda regions saved to {output_path}")
    return nregs, categories[-1]


def lambda_regions_land_ocean_boxes(cfg, output_path, lambdas):
    """Land cells are each their own region; ocean cells are binned into
    1x1 degree boxes, each box its own region. Useful when land-surface
    flux uncertainty needs cell-level resolution but ocean flux
    uncertainty doesn't (much coarser regions there without losing
    anything meaningful).

    Needs cfg.input_files_dynamics_grid_filename and an EXTPAR file with
    FR_LAND on the same grid, as cfg.input_files_extpar_filename.
    """
    with xr.open_dataset(cfg.input_files_dynamics_grid_filename) as dsgrid, \
         xr.open_dataset(cfg.input_files_extpar_filename) as extpar:
        fr_land = extpar["FR_LAND"].squeeze().values
        clon = np.rad2deg(dsgrid["clon"].values)
        clat = np.rad2deg(dsgrid["clat"].values)

    n_cells = fr_land.size
    is_land = fr_land > 0
    box_deg = 1.0

    # Region IDs must be contiguous 1..n_total (used as a 1-based Fortran
    # array index downstream). np.unique's return_inverse gives exactly
    # that once shifted from 0-based.
    keys = np.where(
        is_land[:, None],
        np.c_[np.arange(n_cells) - 10**12,
              np.zeros(n_cells)],  # unique per land cell
        np.c_[np.floor(clon / box_deg),
              np.floor(clat / box_deg)])  # shared per ocean box
    _, reg0 = np.unique(keys, axis=0, return_inverse=True)
    reg = (reg0 + 1).astype(np.int32)

    categories = np.arange(1, len(lambdas) + 1)
    ds_cells = xr.Dataset(
        data_vars={
            "REG": (["cell"], reg),
            "Lambda_indicies": (["cat"], lambdas),
        },
        coords={
            "cell": (["cell"], np.arange(1, n_cells + 1)),
            "cat": (["cat"], categories),
        },
        attrs={"author": "Processing Chain"},
    )
    ds_cells.to_netcdf(output_path,
                       encoding={
                           "REG": {
                               "dtype": "int32"
                           },
                           "cell": {
                               "dtype": "int32"
                           }
                       })
    print(f"Lambda regions saved to {output_path}")
    return int(reg.max()), categories[-1]


def lambda_regions_parent_grid(cfg, output_path, lambdas):
    """Regions = the case's coarser 'parent' grid. Each dynamics-grid cell
    is assigned to its nearest parent-grid cell, so the fine grid's cells
    collapse into however many regions the parent grid has, instead of
    either one region per fine cell or a domain-specific scheme like
    land/ocean binning.

    Needs cfg.input_files_dynamics_grid_filename and
    cfg.input_files_radiation_grid_filename (or any other coarser grid
    file with clon/clat) as the parent grid.
    """
    with xr.open_dataset(cfg.input_files_dynamics_grid_filename) as ds_fine, \
         xr.open_dataset(cfg.input_files_radiation_grid_filename) as ds_parent:
        fine_lon = np.rad2deg(ds_fine["clon"].values)
        fine_lat = np.rad2deg(ds_fine["clat"].values)
        parent_lon = np.rad2deg(ds_parent["clon"].values)
        parent_lat = np.rad2deg(ds_parent["clat"].values)

    n_cells = fine_lon.size
    n_parent = parent_lon.size

    tree = cKDTree(np.column_stack([parent_lon, parent_lat]))
    _, reg0 = tree.query(np.column_stack([fine_lon, fine_lat]))
    reg = (reg0 + 1).astype(np.int32)  # region IDs are 1-based

    categories = np.arange(1, len(lambdas) + 1)
    ds_cells = xr.Dataset(
        data_vars={
            "REG": (["cell"], reg),
            "Lambda_indicies": (["cat"], lambdas),
        },
        coords={
            "cell": (["cell"], np.arange(1, n_cells + 1)),
            "cat": (["cat"], categories),
        },
        attrs={"author": "Processing Chain"},
    )
    ds_cells.to_netcdf(output_path,
                       encoding={
                           "REG": {
                               "dtype": "int32"
                           },
                           "cell": {
                               "dtype": "int32"
                           }
                       })
    print(f"Lambda regions saved to {output_path}")
    return n_parent, categories[-1]


def boundary_regions_angular_quadrant(cfg, output_path, n_bg_ens):
    """The default: split the domain into n_bg_ens angular sectors measured
    from its centroid (via iconsub, needs cfg.cdo_nco_cmd/cdo_nco_cmd_post
    to load cdo/nco/icontools). Simple and fast, but for an elongated or
    irregular domain the sectors can come out badly lopsided, since
    "nearest boundary segment" isn't the same thing as "narrowest angular
    slice" once the domain isn't roughly circular -- see
    boundary_regions_compass_walk for an alternative that follows the
    domain's actual shape instead.

    Needs cfg.input_files_dynamics_grid_filename, cfg.cdo_nco_cmd,
    cfg.cdo_nco_cmd_post, cfg.user_name, cfg.user_mail.
    """
    grid_filename = cfg.input_files_dynamics_grid_filename
    workdir = output_path.parent / "iconsub_work"
    workdir.mkdir(parents=True, exist_ok=True)

    cmd = f"""
{cfg.cdo_nco_cmd}
cat > NAMELIST_ICONSUB << EOF_1
&iconsub_nml
    grid_filename = '{grid_filename}',
    output_type = 4,
    lwrite_grid = .TRUE.,
/
&subarea_nml
    ORDER = "outgrid",
    grf_info_file = '{grid_filename}',
    min_refin_c_ctrl = 1,
    max_refin_c_ctrl = 120
/
EOF_1

iconsub --nml NAMELIST_ICONSUB
{cfg.cdo_nco_cmd_post}
    """
    subprocess.check_output(cmd, shell=True, cwd=workdir)

    ds_grid = xr.open_dataset(workdir / 'outgrid.grid.nc')
    clon, clat = np.rad2deg(ds_grid['clon']), np.rad2deg(ds_grid['clat'])

    # Compute the central reference point
    mid_lon, mid_lat = np.nanquantile(clon, 0.5), np.nanquantile(clat, 0.5)

    # Center coordinates relative to the midpoint
    clon_cent, clat_cent = clon - mid_lon, clat - mid_lat

    # Compute angles of all points relative to the center
    angles = np.arctan2(clat_cent, clon_cent)  # Range: [-π, π]

    # Set number of regions
    sector_size = (2 * np.pi) / n_bg_ens  # Each sector covers an angle range

    # Assign each point to a region (0 to N-1)
    region_indices = (angles // sector_size).astype(int)

    # One-hot encode the region assignments
    boundary_regions = np.zeros((len(clon), n_bg_ens), dtype=np.int32)
    boundary_regions[np.arange(len(clon)), region_indices] = 1

    attrs = {'author': cfg.user_name}
    if cfg.user_mail:
        attrs['email'] = cfg.user_mail
    ds_boundary = xr.Dataset(data_vars={
        'boundaryregion': (['cell', 'reg'], boundary_regions),
        'global_cell_idx': (['cell'], np.arange(len(clon)))
    },
                             coords={
                                 'cell': (['cell'], np.arange(len(clon))),
                                 'reg': (['reg'], np.arange(n_bg_ens))
                             },
                             attrs=attrs)
    try:
        ds_boundary.to_netcdf(output_path)
    except:
        print("File currently open. Please close the file and try again.")
    print(f"Boundary regions saved to {output_path}")


def boundary_regions_compass_walk(cfg, output_path, n_bg_ens):
    """Walk the domain's actual boundary (via grid neighbor connectivity)
    into an ordered loop, split it into 8 compass-direction segments
    (NNE/ENE/ESE/SSE/SSW/WSW/WNW/NNW) anchored on the bounding box's
    N/NE/E/SE/S/SW/W/NW points, then assign every interior cell to its
    nearest boundary segment. Unlike boundary_regions_angular_quadrant,
    this follows the domain's actual shape.

    Only supports n_bg_ens == 8 (one region per compass direction) -- the
    segment-splitting logic below is written around exactly 8 named
    anchors.

    Needs cfg.input_files_dynamics_grid_filename (with a
    neighbor_cell_index variable), cfg.user_name, cfg.user_mail.

    Note on indexing: this follows boundary_regions_angular_quadrant's
    convention of 0-based cell/reg coordinates in the output NetCDF (the
    lambda-regions strategies above are 1-based instead, since their REG
    variable is used directly as a 1-based Fortran index downstream --
    boundaryregion here is a one-hot matrix, not an index value, so the
    coordinate numbering itself doesn't carry that same constraint).
    """
    if n_bg_ens != 8:
        raise NotImplementedError(
            "This compass-segment boundary-region scheme only supports "
            "n_bg_ens == 8 (one region per compass direction)")

    with xr.open_dataset(cfg.input_files_dynamics_grid_filename) as ds:
        clon = np.rad2deg(ds['clon'].values)
        clat = np.rad2deg(ds['clat'].values)
        ncell = clon.size
        raw_nbr = ds['neighbor_cell_index'].values

    if raw_nbr.ndim != 2:
        raise RuntimeError(
            "Unexpected neighbor array shape: expected 2D array.")
    if raw_nbr.shape[0] == ncell:
        nbr = raw_nbr.copy()
    elif raw_nbr.shape[1] == ncell:
        nbr = raw_nbr.T.copy()
    else:
        raise RuntimeError(
            f"Can't interpret neighbor array shape {raw_nbr.shape} for ncell={ncell}"
        )

    # ICON grid files use either 0 or 1-based neighbor indices with 0/negative
    # meaning "no neighbor" (domain edge) -- normalize both to 0-based with -1
    # for missing.
    nbr = nbr.astype(int)
    if (nbr <= 0).any() and (nbr.min() == 0):
        nbr = np.where(nbr <= 0, -1, nbr - 1)
    elif nbr.min() == 1:
        nbr = nbr - 1

    # Boundary cells: at least one missing neighbor
    boundary_mask = np.any(nbr < 0, axis=1)
    boundary_cell_idx = np.where(boundary_mask)[0]
    b_clon = clon[boundary_cell_idx]
    b_clat = clat[boundary_cell_idx]
    boundary_coords = np.column_stack([b_clon, b_clat])

    # Walk the boundary into an ordered loop via adjacency restricted to
    # boundary cells, starting from the westernmost one.
    boundary_set = set(boundary_cell_idx.tolist())
    adj = {ci: [] for ci in boundary_cell_idx}
    for ci in boundary_cell_idx:
        for nb in nbr[ci]:
            if nb >= 0 and nb in boundary_set:
                adj[ci].append(int(nb))

    start = int(boundary_cell_idx[np.argmin(b_clon)])
    loop = []
    visited = set()
    curr, prev = start, None
    max_steps = len(boundary_cell_idx) + 50
    for _ in range(max_steps):
        loop.append(curr)
        visited.add(curr)
        cand = [nb for nb in adj[curr] if nb != prev]
        unvisited = [nb for nb in cand if nb not in visited]
        if unvisited:
            next_cell = unvisited[0]
        elif cand:
            next_cell = cand[0]
        else:
            break
        prev, curr = curr, next_cell
        if curr == start:
            break
    boundary_idx = np.array(loop, dtype=int)

    # A degenerate/broken walk (implausibly short) falls back to a simple
    # polar sort around the boundary's centroid.
    if len(boundary_idx) < max(10, int(0.02 * len(boundary_cell_idx))):
        centroid = np.array(
            [boundary_coords[:, 0].mean(), boundary_coords[:, 1].mean()])
        ang = np.arctan2(boundary_coords[:, 1] - centroid[1],
                         boundary_coords[:, 0] - centroid[0])
        boundary_idx = boundary_cell_idx[np.argsort(ang)]
    boundary_coords = np.column_stack([clon[boundary_idx], clat[boundary_idx]])

    # Bounding-box anchors: N, NE, E, SE, S, SW, W, NW
    minx, maxx = clon.min(), clon.max()
    miny, maxy = clat.min(), clat.max()
    anchors = np.array([
        [0.5 * (minx + maxx), maxy],
        [maxx, maxy],
        [maxx, 0.5 * (miny + maxy)],
        [maxx, miny],
        [0.5 * (minx + maxx), miny],
        [minx, miny],
        [minx, 0.5 * (miny + maxy)],
        [minx, maxy],
    ])

    l_loop = boundary_coords.shape[0]
    anchor_idx_in_loop = [
        int(
            np.argmin(
                np.hypot(boundary_coords[:, 0] - ax,
                         boundary_coords[:, 1] - ay))) for ax, ay in anchors
    ]
    anchor_idx_sorted = np.sort(np.array(anchor_idx_in_loop, dtype=int))

    # Split the ordered loop into the 8 named segments between consecutive
    # (sorted) anchors, then map every boundary cell to its segment index.
    boundary_segment_of_cell = -1 * np.ones(ncell, dtype=int)
    for i in range(8):
        start_i = anchor_idx_sorted[i]
        end_i = anchor_idx_sorted[(i + 1) % 8]
        if end_i > start_i:
            seg_loop_pos = np.arange(start_i, end_i + 1)
        else:
            seg_loop_pos = np.concatenate(
                [np.arange(start_i, l_loop),
                 np.arange(0, end_i + 1)])
        boundary_segment_of_cell[boundary_idx[seg_loop_pos]] = i

    # Assign every interior cell to its nearest boundary segment.
    interior_cells = np.where(boundary_segment_of_cell == -1)[0]
    cell_xy = np.column_stack([clon, clat])
    all_boundary_cells = np.where(boundary_segment_of_cell != -1)[0]
    bxy = cell_xy[all_boundary_cells]
    bseg = boundary_segment_of_cell[all_boundary_cells]

    cell_region = boundary_segment_of_cell.copy()
    batch = 5000
    for batch_start in range(0, len(interior_cells), batch):
        idx = interior_cells[batch_start:batch_start + batch]
        xy = cell_xy[idx]
        dx = xy[:, None, 0] - bxy[None, :, 0]
        dy = xy[:, None, 1] - bxy[None, :, 1]
        nearest = np.argmin(dx * dx + dy * dy, axis=1)
        cell_region[idx] = bseg[nearest]

    attrs = {"author": cfg.user_name}
    if cfg.user_mail:
        attrs["email"] = cfg.user_mail
    ds_boundary = xr.Dataset(
        data_vars={
            "boundaryregion":
            (["cell", "reg"], np.eye(8, dtype=np.int32)[cell_region]),
            "global_cell_idx": (["cell"], np.arange(ncell)),
        },
        coords={
            "cell": (["cell"], np.arange(ncell)),
            "reg": (["reg"], np.arange(8)),
        },
        attrs=attrs,
    )
    ds_boundary.to_netcdf(output_path)
    print(f"Boundary regions saved to {output_path}")


LAMBDA_REGION_STRATEGIES = {
    'per_cell': lambda_regions_per_cell,
    'land_ocean_boxes': lambda_regions_land_ocean_boxes,
    'parent_grid': lambda_regions_parent_grid,
}

BOUNDARY_REGION_STRATEGIES = {
    'angular_quadrant': boundary_regions_angular_quadrant,
    'compass_walk': boundary_regions_compass_walk,
}
