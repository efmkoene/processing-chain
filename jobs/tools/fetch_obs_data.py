import os
import re
import logging
from pathlib import Path
from datetime import datetime, timedelta
from time import sleep
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import xarray as xr
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
from shapely.geometry import Polygon
from shapely import contains_xy
from icoscp.dobj import Dobj
from icoscp.sparql.runsparql import RunSparql

from . import iter_hours


def fetch_ICOS_data(query_type='any',
                    start_date='01-01-2022',
                    end_date='31-12-2022',
                    save_path='',
                    species=['co', 'co2', 'ch4'],
                    collection_url=None):
    """
    Runs a SPARQL query against the ICOS-CP portal and downloads the
    matching data objects as one NetCDF file per station.

    Parameters
    ----------
    query_type : str
        One of 'release', 'growing', 'any' -- selects between the different
        file products at the ICOS-CP.
    start_date, end_date : str
        Dates in dd-mm-yyyy format.
    save_path : str
        Directory to save the downloaded per-station files to.
    species : list of str
        Subset of ['co', 'co2', 'ch4'].
    collection_url : str, optional
        If given, download every member of this ICOS-CP collection instead
        of running the spec-based (atc{species}Product) search above.
        Needed for non-ICOS-network partner stations that only exist in
        curated ObsPack collections (spec ObspackTimeSerieResult), e.g. the
        European CO2 ObsPack compilation at
        https://meta.icos-cp.eu/collections/-UpVwVSamTSg-HEDhh5oz-0d
    """
    if collection_url is not None:
        query = '''
        prefix cpmeta: <http://meta.icos-cp.eu/ontologies/cpmeta/>
        prefix dcterms: <http://purl.org/dc/terms/>
        select ?dobj ?fileName where {{
            <{0}> dcterms:hasPart ?dobj .
            ?dobj cpmeta:hasName ?fileName .
        }}
        order by ?fileName
        '''.format(collection_url)
    else:
        # --- Build up an SQL query for the different species
        qd = ""
        for specie in species:
            qd += f" <http://meta.icos-cp.eu/resources/cpmeta/atc{specie.capitalize()}"
            if query_type == 'release':
                qd += "L2DataObject>"
            elif query_type == 'growing':
                qd += "NrtGrowingDataObject>"
            elif query_type == 'any':
                qd += "Product>"

        query = '''
        prefix cpmeta: <http://meta.icos-cp.eu/ontologies/cpmeta/>
        prefix prov: <http://www.w3.org/ns/prov#>
        prefix xsd: <http://www.w3.org/2001/XMLSchema#>
        select ?dobj ?hasNextVersion ?spec ?fileName ?size ?submTime ?timeStart ?timeEnd
        where {{
            VALUES ?spec {{{0}}}
            ?dobj cpmeta:hasObjectSpec ?spec .
            BIND(EXISTS{{[] cpmeta:isNextVersionOf ?dobj}} AS ?hasNextVersion)
            ?dobj cpmeta:hasSizeInBytes ?size .
        ?dobj cpmeta:hasName ?fileName .
        ?dobj cpmeta:wasSubmittedBy/prov:endedAtTime ?submTime .
        ?dobj cpmeta:hasStartTime | (cpmeta:wasAcquiredBy / prov:startedAtTime) ?timeStart .
        ?dobj cpmeta:hasEndTime | (cpmeta:wasAcquiredBy / prov:endedAtTime) ?timeEnd .
            FILTER NOT EXISTS {{[] cpmeta:isNextVersionOf ?dobj}}
        FILTER( !(?timeStart > '{1}T23:00:00.000Z'^^xsd:dateTime || ?timeEnd < '2017-12-31T23:00:00.000Z'^^xsd:dateTime) )

        }}
        order by desc(?submTime)
        '''.format(qd, (datetime.strptime(start_date, '%d-%m-%Y').date() -
                        timedelta(days=1)).strftime('%Y-%m-%d'),
                   (datetime.strptime(end_date,
                                      '%d-%m-%Y').date()).strftime('%Y-%m-%d'))

    # --- Run the SQL query
    result = RunSparql(query, 'pandas')
    result.run()
    result.data()

    # --- Loop over the different stations (see https://icos-carbon-portal.github.io/pylib/ for more details)
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    specie = set(species)  # placeholder until the per-object line below narrows it;
                           # needed unconditionally since the collection_url path
                           # never runs the species-labeled spec-query loop above
    for d in result.data()['dobj']:
        logging.info(f"Processing data object: {specie} for object {d}")
        finished = False
        while not finished:
            try:
                outfn = os.path.join(
                    save_path, 'ICOS_obs_' + str(specie)[2:-2] + '_' + query_type +
                    '_' + str(Dobj(d).station['id']) + '_' +
                    str(Dobj(d).meta['specificInfo']['acquisition']['samplingHeight'])
                    + '_' + start_date + '_' + end_date + '.nc')

                # Skip if filename exists (checked before the expensive .data pull,
                # so a re-run over already-cached stations doesn't re-fetch them)
                if os.path.isfile(outfn):
                    finished = True
                    continue
                obj = Dobj(d).data

                lon = Dobj(d).lon
                lat = Dobj(d).lat
                variables = Dobj(d).variables.to_numpy()
                Names = Dobj(d).colNames
                specie = set(Names) - set(Names).difference(species)
                meta = np.squeeze(
                    [x for x in variables if set(species) - set(x) != set(species)])
                ds = xr.Dataset.from_dataframe(obj)  # This contains the data...
                # --- Cleanup of the dataframe...
                ds = ds.set_index(index='TIMESTAMP')
                ds = ds.sortby(ds.index)
                ds = ds.drop_duplicates(dim="index")
                # --- Subset to the timeframe of interest (this has no reason to fail, so you'll have to check these cases manually....)
                ds = ds.sel(index=slice(
                    datetime.strptime(start_date, '%d-%m-%Y').date().strftime(
                        '%Y-%m-%d'),
                    datetime.strptime(end_date, '%d-%m-%Y').date().strftime(
                        '%Y-%m-%d')))

                ds = ds.rename({'index': 'time'})
                # --- Write out further attributes
                ds.attrs['Description'] = meta[2]
                ds.attrs['Units'] = meta[1]
                ds.attrs['Station'] = Dobj(d).station['id']
                ds.attrs['Full name of the station'] = Dobj(d).station['org']['name']
                ds.attrs['Elevation above sea level'] = Dobj(d).alt
                ds.attrs['Sampling height over ground'] = Dobj(
                    d).meta['specificInfo']['acquisition']['samplingHeight']
                ds.attrs['Sampling height over sea level'] = float(
                    Dobj(d).meta['specificInfo']['acquisition']
                    ['samplingHeight']) + float(Dobj(d).alt)
                ds.attrs['Longitude'] = Dobj(d).lon
                ds.attrs['Latitude'] = Dobj(d).lat
                ds.attrs['Name of the tracer'] = meta[0]
                var_name = next(iter(specie))
                if not np.any(np.isfinite(ds[var_name].values)):
                    # Station genuinely has no valid data in this window
                    # (e.g. before it came online) -- process_ICOS_data
                    # would exclude it anyway (any(np.isfinite(cnc)) check),
                    # but writing a placeholder file here under the same
                    # name a *real* future window's data would use blocks
                    # that later fetch via the if-exists skip below.
                    logging.info(
                        f"No valid {var_name} data for {os.path.basename(outfn)} "
                        f"in this window, not writing a placeholder file")
                else:
                    ds.to_netcdf(outfn)
                finished = True
            except Exception:
                logging.info('Waiting for server to respond...')
                sleep(5)


def process_ICOS_data(ICOS_obs_folder,
                      start_date='01-01-2022',
                      end_date='31-12-2022',
                      output_folder='~/',
                      icon_poly=None):
    """Package the downloaded per-station ICOS data into a single file.

    Parameters
    ----------
    ICOS_obs_folder : str
        Directory the per-station files from `fetch_ICOS_data` were saved to.
    start_date, end_date : datetime
    output_folder : str
    icon_poly : shapely Polygon/MultiPolygon, optional
        Domain to test each station against with contains_xy() (e.g. from
        `load_icon_polygon_erode_buffer()`). Falls back to the lon_lims/
        lat_lims bounding box below if not given.
    """
    output_filename = Path(
        output_folder
    ) / f"Extracted_{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}_alldates_masl.nc"
    if os.path.isfile(output_filename):
        return

    # Fallback bounding box, only used if icon_poly isn't given. A box is
    # cruder than icon_poly's actual shape (over-includes stations in the
    # box's corners that are outside the true, irregular domain), so pass
    # icon_poly when possible.
    lon_lims = [-10.76, 22.98]
    lat_lims = [36.69, 60.35]

    # Utility for converting units to PPMv
    toppm_dict = {'nmol mol-1': 1e-9 * 1e6, 'µmol mol-1': 1e-6 * 1e6}

    # Gather chosen dates
    delta = end_date - start_date
    chosen_dates = [
        np.datetime64((start_date + timedelta(
            days=i, hours=h)).strftime('%Y-%m-%dT%H:%M:%S.000000000'))
        for i in range(delta.days + 1) for h in range(24)
    ]
    number_of_hourly_measurements = len(chosen_dates)
    logging.info(
        f'A total of {number_of_hourly_measurements} hours are possible')

    # Gather files
    logging.info(
        f"Looking in folder {ICOS_obs_folder} for ICOS observation files with glob *{start_date.strftime('%d-%m-%Y')}_{end_date.strftime('%d-%m-%Y')}.nc"
    )
    files = list(
        Path(ICOS_obs_folder).glob(
            f"*{start_date.strftime('%d-%m-%Y')}_{end_date.strftime('%d-%m-%Y')}.nc"
        ))
    number_of_stations = len(files)
    logging.info(f'Will package data from {number_of_stations} files, {files}')

    if number_of_stations == 0:
        logging.info('No ICOS station files found for this window, nothing to package')
        return

    # Prepare
    obs_cnc_matrix = np.zeros(
        (number_of_stations, number_of_hourly_measurements), dtype=np.float64)
    obs_dates_matrix = np.zeros(
        (number_of_stations, number_of_hourly_measurements),
        dtype=np.dtype('datetime64[ns]'))
    obs_std_matrix = np.zeros(
        (number_of_stations, number_of_hourly_measurements), dtype=np.float64)

    # Set-up a function that can be called in parallel
    def extract_obs_column(file):
        logging.info(f'Opened file {file}')
        try:
            # Open dataset and extract metadata
            ds = xr.open_dataset(file)
            name = f"{ds.attrs['Full name of the station']}_{file.name.split('_')[-3][:-2]}"
            id_st = ds.attrs['Station']
            units = ds.attrs['Units']
            masl = ds.attrs['Elevation above sea level']
            diff = (ds.time.values[1] - ds.time.values[0]
                    ) / 3600000000000  # Time difference in hours

            if diff != 1:
                logging.info(
                    f'Observation data at station {name} is not hourly averaged ({diff} hours)'
                )

            # Filter dataset to the desired time range. Slicing with bare
            # datetimes (rather than date strings) is an exact-timestamp
            # bound, not an inclusive whole-day one -- it silently drops
            # all but the midnight sample of end_date. Use date strings, as
            # fetch_ICOS_data already does for the same reason above.
            ds_filtered = ds.sel(time=slice(
                start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')))

            # Align `chosen_dates` with `ds_filtered.time`
            ds_aligned = ds_filtered.reindex(time=chosen_dates,
                                             method='nearest',
                                             tolerance='1h')

            # Update observation arrays
            obs_dates1 = ds_aligned.time.values
            obs_std1 = ds_aligned.Stdev.values * toppm_dict[units]
            obs_cnc1 = ds_aligned["co2"].values * toppm_dict[units]
            lons, lats = ds.attrs['Longitude'], ds.attrs['Latitude']

        except Exception as e:
            logging.info(f"Error processing file {file}: {e}")
            obs_cnc1 = np.full(number_of_hourly_measurements,
                               np.nan,
                               dtype=np.float64)
            obs_dates1 = np.full(number_of_hourly_measurements,
                                 np.datetime64("NaT"),
                                 dtype="datetime64[ns]")
            obs_std1 = np.full(number_of_hourly_measurements,
                               np.nan,
                               dtype=np.float64)
            name, id_st, masl, lons, lats = 'nan', 0, -999, np.nan, np.nan

        return name, obs_std1, obs_cnc1, obs_dates1, lons, lats, id_st, masl

    # Process all data concurrently
    with ThreadPoolExecutor(max_workers=1) as executor:
        results = list(executor.map(extract_obs_column, files))
    M = list(zip(*results))

    station_names = np.array(M[0])
    obs_cnc = np.array(M[2])
    obs_std = np.array(M[1])
    obs_times = np.array(M[3])
    obs_lons = np.array(M[4])
    obs_lats = np.array(M[5])
    obs_ids = np.array(M[6])
    obs_masl = np.array(M[7])

    # Initialize mask and removal list
    stations_to_keep = []
    mask_true = np.full_like(obs_cnc_matrix[0], True)

    # Filter and populate matrices
    for ix, (lon, lat, cnc, std, times) in enumerate(
            zip(obs_lons, obs_lats, obs_cnc, obs_std, obs_times)):
        if icon_poly is not None:
            in_domain = bool(contains_xy(icon_poly, lon, lat))
        else:
            in_domain = (lon_lims[0] < lon < lon_lims[-1]) and (lat_lims[0] < lat < lat_lims[-1])
        if any(np.isfinite(cnc)) and in_domain:
            np.place(obs_cnc_matrix[ix], mask_true, cnc)
            np.place(obs_std_matrix[ix], mask_true, std)
            np.place(obs_dates_matrix[ix], mask_true, times)
            stations_to_keep.append(ix)

    # Convert keep list to numpy index array for slicing
    stations_to_keep = np.array(stations_to_keep)

    # Filter matrices and metadata
    obs_cnc_matrix = obs_cnc_matrix[stations_to_keep]
    obs_std_matrix = obs_std_matrix[stations_to_keep]
    obs_dates_matrix = obs_dates_matrix[stations_to_keep]
    station_names = station_names[stations_to_keep]
    obs_lons = obs_lons[stations_to_keep]
    obs_lats = obs_lats[stations_to_keep]
    obs_ids = obs_ids[stations_to_keep]
    obs_masl = obs_masl[stations_to_keep]
    station_idcs = np.arange(len(station_names))

    # Define data variables and attributes for xarray dataset
    data_vars = {
        "Concentration": (["station", "time"], obs_cnc_matrix, {
            "units": "ppm",
            "long_name": "CO2_concentration"
        }),
        "Std": (["station", "time"], obs_std_matrix, {
            "units": "ppm",
            "long_name": "CO2_concentrations_std"
        }),
        "Stations_names": (["station"], station_names, {
            "units": "-",
            "long_name": "Stations_names"
        }),
        "Stations_ids": (["station"], obs_ids, {
            "units": "-",
            "long_name": "Stations_names"
        }),
        "Stations_masl": (["station"], obs_masl, {
            "units": "-",
            "long_name": "Elevation_heights_above_sl"
        }),
        "Lon": (["station"], obs_lons, {
            "units": "degrees",
            "long_name": "Longitude"
        }),
        "Lat": (["station"], obs_lats, {
            "units": "degrees",
            "long_name": "Latitude"
        }),
        "Dates": (["station", "time"], obs_dates_matrix, {
            "long_name": "Dates"
        }),
    }

    # Define coordinates
    coords = {"station": (["station"], station_idcs)}
    attrs = {
        'creation_date': str(datetime.now()),
        'author': 'Processing Chain'
    }

    # Create xarray dataset
    ds_extracted_obs_matrix = xr.Dataset(data_vars=data_vars,
                                         coords=coords,
                                         attrs=attrs)

    # Save dataset to file
    ds_extracted_obs_matrix.to_netcdf(output_filename)

    logging.info(
        f"Finished extraction and stored obs_matrix for {len(obs_lons)} stations "
        f"(from {number_of_stations} available ICOS stations), which were operating "
        f"during the given period and are located inside the model domain, in the file: {output_filename}"
    )


def list_nc4_urls(index_url, DATE_WINDOW, session):
    """Scrape index_url and return absolute URLs ending with .nc4 (not .nc4.xml)."""
    logging.info(f"Fetching directory listing: {index_url}")
    r = session.get(index_url)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, "html.parser")
    links = [a.get("href") for a in soup.find_all("a", href=True)]
    urls = []
    NC4_REGEX = re.compile(r".*\.nc4$", re.IGNORECASE)
    for href in links:
        if NC4_REGEX.match(href) and not href.endswith(".nc4.xml"):
            # resolve relative URLs
            if href.startswith("http"):
                urls.append(href)
            else:
                urls.append(requests.compat.urljoin(index_url, href))
    # optional filename date filter using tokens like yyyymmdd or yymmdd in filename
    if DATE_WINDOW:
        sdt = datetime.strptime(DATE_WINDOW[0], "%Y-%m-%d").date()
        edt = datetime.strptime(DATE_WINDOW[1], "%Y-%m-%d").date()

        def date_in_name(u):
            m = re.search(r"(\d{8}|\d{6})", u)
            if not m:
                return True
            token = m.group(1)
            try:
                if len(token) == 8:
                    d = datetime.strptime(token, "%Y%m%d").date()
                else:
                    d = datetime.strptime(token, "%y%m%d").date()
            except Exception:
                return True
            return sdt <= d <= edt

        urls = [u for u in urls if date_in_name(u)]
    urls.sort()
    logging.info(f"Found {len(urls)} .nc4 files (after optional filtering).")
    return urls


def download_file(url, out_folder: Path, session=None):
    """Download URL into out_folder, streaming with progress bar. Uses session (honors .netrc auth)."""
    session = session or requests.Session()
    local = out_folder / Path(url).name
    if local.exists():
        logging.info(f"Skipping download (exists): {local.name}")
        return local
    with session.get(url, stream=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(local, "wb") as f, tqdm(total=total, unit="B", unit_scale=True, desc=local.name) as pbar:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
    return local


def load_icon_polygon_erode_buffer(icon_nc_path, n_layers=42):
    """
    Load ICON unstructured grid and erode it by n_layers along the true
    geometric boundary determined from unique edges.

    Args:
        icon_nc_path (str or Path): Path to icon_domain.nc
        n_layers (int): Number of layers to erode

    Returns:
        Shapely Polygon/MultiPolygon of the eroded domain
    """
    ds = xr.open_dataset(icon_nc_path)

    # --- Build cell polygons ---
    clon_v = ds["clon_vertices"].values
    clat_v = ds["clat_vertices"].values

    if np.nanmax(np.abs(clon_v)) < 4:
        clon_v = np.degrees(clon_v)
        clat_v = np.degrees(clat_v)

    n_cells = clon_v.shape[0]
    polys = []
    for i in range(n_cells):
        coords = np.column_stack((clon_v[i, :], clat_v[i, :]))
        coords = coords[~np.isnan(coords).any(axis=1)]
        if coords.shape[0] >= 3:
            polys.append(Polygon(coords))
        else:
            polys.append(None)

    ds.close()

    # --- Iterative erosion by geometric boundary ---
    active = set(range(n_cells))
    for layer in range(n_layers):
        # Recompute edges among active cells
        active_edge_count = defaultdict(int)
        edge_to_cells = defaultdict(list)

        for idx in active:
            poly = polys[idx]
            if poly is None:
                continue
            coords = np.array(poly.exterior.coords)
            for i in range(len(coords) - 1):
                edge = tuple(sorted([tuple(coords[i]), tuple(coords[i + 1])]))
                active_edge_count[edge] += 1
                edge_to_cells[edge].append(idx)

        # Identify boundary cells: any cell that has at least one edge appearing only once
        boundary = set()
        for edge, count in active_edge_count.items():
            if count == 1:
                boundary.add(edge_to_cells[edge][0])

        if not boundary:
            break

        active -= boundary
        logging.info(f"Layer {layer + 1}: {len(boundary)} cells removed, {len(active)} remaining")

    # --- Union remaining polygons ---
    # shapely.ops.unary_union() routes through shapely's vectorized
    # union_all ufunc, which breaks under some numpy/shapely combinations
    # ("ufunc 'create_collection' not supported for the input types").
    # Binary Polygon.union() doesn't go through that ufunc, so a
    # tree-reduction (pairwise, halving each round) union avoids it while
    # staying fast.
    eroded_polys = [polys[i] for i in active if polys[i] is not None]

    def _tree_union(geoms):
        geoms = list(geoms)
        while len(geoms) > 1:
            nxt = [geoms[i].union(geoms[i + 1]) for i in range(0, len(geoms) - 1, 2)]
            if len(geoms) % 2:
                nxt.append(geoms[-1])
            geoms = nxt
        return geoms[0]

    domain = _tree_union(eroded_polys)

    logging.info(f"Eroded ICON polygon built ({domain.geom_type}, area={domain.area:.4f} deg^2)")
    return domain


def normalize_longitudes(lon):
    """Convert lon array to range [-180, 180]"""
    return ((lon + 180) % 360) - 180


def mask_dataset_to_polygon(ds, poly):
    # Find lat/lon
    lat_name = next((v for v in ds if "lat" in v.lower()), None)
    lon_name = next((v for v in ds if "lon" in v.lower()), None)
    if lat_name is None or lon_name is None:
        raise ValueError("Could not identify lat/lon variables")

    lat = ds[lat_name].values
    lon = normalize_longitudes(ds[lon_name].values)

    # Determine which soundings are inside the polygon
    inside = contains_xy(poly, lon.ravel(), lat.ravel())
    mask = inside.reshape(lat.shape[0])

    return ds.isel(sounding_id=mask)


def fetch_OCO2_data(DATE_WINDOW=("2017-12-28", "2018-01-15"),
                    OUT_DIR=Path("./oco2_downloads"),
                    ICON_GRID_PATH="icon_domain.nc",
                    n_layers=12):
    """
    Fetch OCO-2 L2 Lite data from NASA GES DISC, mask to the ICON domain,
    and save the masked files.

    Parameters
    ----------
    DATE_WINDOW : tuple of str
        (start_date, end_date) in "YYYY-MM-DD" format.
    OUT_DIR : Path
        Output directory for downloaded and masked files.
    ICON_GRID_PATH : str
        Path to ICON grid netCDF file.
    n_layers : int
        Number of layers to erode the ICON domain by.
    """
    logging.info(
        f"Fetch OCO2 called with DATE_WINDOW={DATE_WINDOW}, OUT_DIR={OUT_DIR}, ICON_GRID_PATH={ICON_GRID_PATH}, n_layers={n_layers}"
    )

    session = requests.Session()

    OUT_DIR = Path(OUT_DIR)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    start_year = int(DATE_WINDOW[0][:4])
    end_year = int(DATE_WINDOW[1][:4])
    years = list(range(start_year, end_year + 1))

    BASE_URL_TEMPLATE = "https://oco2.gesdisc.eosdis.nasa.gov/data/OCO2_DATA/OCO2_L2_Lite_FP.11.2r/{year}/"
    base_dir_listings = [BASE_URL_TEMPLATE.format(year=y) for y in years]

    all_urls = []
    for base_dir in base_dir_listings:
        urls = list_nc4_urls(base_dir, DATE_WINDOW, session)
        all_urls.extend(urls)
    if not all_urls:
        logging.info("No .nc4 files found. Exiting.")
        return

    # Final product is the masked file, so that's what determines whether a
    # URL needs anything done at all -- skips both the (potentially large)
    # raw download and the masking/write for URLs already processed by a
    # previous run.
    def masked_path(url):
        return OUT_DIR / (Path(url).stem + "_masked.nc4")

    urls_to_process = [u for u in all_urls if not masked_path(u).exists()]
    masked_files = [masked_path(u) for u in all_urls if masked_path(u).exists()]
    if not urls_to_process:
        logging.info(
            f"All {len(all_urls)} OCO-2 files already downloaded and masked, nothing to do"
        )
        return

    logging.info("Loading ICON domain and building domain polygon (this may take some seconds)...")
    icon_poly = load_icon_polygon_erode_buffer(ICON_GRID_PATH, n_layers=n_layers)
    logging.info("ICON polygon built.")

    skipped = []
    for url in urls_to_process:
        try:
            local = download_file(url, OUT_DIR, session)
        except Exception as e:
            logging.info(f"Download failed for {url}: {e}")
            continue
        # open, mask, and save
        ds = xr.open_dataset(str(local))
        try:
            ds_masked = mask_dataset_to_polygon(ds, icon_poly)
        except Exception as e:
            logging.info(f"Masking failed for {local.name}: {e}")
            ds.close()
            continue
        out_masked = masked_path(url)
        logging.info(f"Saving masked file: {out_masked.name}")
        enc = {v: {"zlib": True, "complevel": 4} for v in ds_masked.data_vars}
        ds_masked.to_netcdf(str(out_masked), format="NETCDF4", encoding=enc)
        ds.close()
        masked_files.append(out_masked)

    logging.info(f"Processing complete. Masked: {len(masked_files)} files. Skipped: {len(skipped)}.")


def process_OCO2_data(OCO2_obs_folder,
                      ICON_grid_file,
                      start_date='01-01-2022',
                      end_date='31-12-2022',
                      output_folder='~/'):
    """Package the downloaded, masked OCO-2 data into one file per day.

    Parameters
    ----------
    OCO2_obs_folder : str
        Directory the masked files from `fetch_OCO2_data` were saved to.
    ICON_grid_file : str
    start_date, end_date : datetime
    output_folder : str
    """
    output_folder = Path(output_folder)

    def _empty_dataset(retrieval_id):
        return xr.Dataset(
            {
                "latitude": (["soundings"], np.array([], dtype=np.float32)),
                "longitude": (["soundings"], np.array([], dtype=np.float32)),
                "date": (["soundings", "epoch_dimension"
                        ], np.empty((0, 7), dtype=np.float32)),
                "obs": (["soundings"], np.array([], dtype=np.float32)),
                "quality_flag": (["soundings"], np.array([], dtype=np.int32)),
                "averaging_kernel": (["soundings", "layers"
                                    ], np.empty((0, 20), dtype=np.float32)),
                "pressure_levels": (["soundings", "layers"
                                    ], np.empty((0, 20), dtype=np.float32)),
                "pressure_weighting_function":
                (["soundings", "layers"], np.empty((0, 20), dtype=np.float32)),
                "prior_profile": (["soundings", "layers"
                                  ], np.empty((0, 20), dtype=np.float32)),
                "prior": (["soundings"], np.array([], dtype=np.float32)),
                "uncertainty": (["soundings"], np.array([], dtype=np.float32)),
                "surface_pressure": (["soundings"], np.array([], dtype=np.float32)),
            },
            coords={
                "soundings": np.array([], dtype=np.int32),
                "layers": np.arange(20),
                "epoch_dimension": np.arange(7),
            },
            attrs={
                'creation_date': str(datetime.now()),
                'author': 'Processing Chain',
                'level_def': 'pressure_boundaries',
                'retrieval_id': retrieval_id,
            },
        )

    for day in iter_hours(start_date, end_date, 24):
        output_file = output_folder / f"OCO2_{day.strftime('%Y%m%d')}_ctdas.nc"
        if output_file.exists():
            logging.info(
                f"OCO2 output already exists for {day.strftime('%Y-%m-%d')}, skipping: {output_file}"
            )
            continue

        logging.info(
            f"Looking in folder {OCO2_obs_folder} for OCO-2 observation files with glob oco2_LtCO2_{day.strftime('%y%m%d')}_*masked.nc4"
        )
        file = list(
            Path(OCO2_obs_folder).glob(
                f"oco2_LtCO2_{day.strftime('%y%m%d')}_*masked.nc4"))
        logging.info(f'Found file(s): {file}')

        if not file:
            logging.info(f'No OCO-2 files found for date {day.strftime("%Y-%m-%d")}')
            _empty_dataset('unknown').to_netcdf(output_file)
            continue
        elif len(file) > 1:
            raise IndexError("Error, more OCO-2 files exist than expected. Review.")
        else:
            logging.info(f'Will open data from {file}')

        # Open file
        s5p_data = xr.open_dataset(file[0])

        # Limit to extent of ICON grid
        ICON_grid = xr.open_dataset(ICON_grid_file)
        offset = 0  # degrees offset to ensure no data is beyond the grid bounds
        try:
            s5p_data = s5p_data.where(
                (s5p_data.longitude
                 >= np.rad2deg(ICON_grid.clon.min().values) + offset) &
                (s5p_data.longitude
                 <= np.rad2deg(ICON_grid.clon.max().values) - offset) &
                (s5p_data.latitude
                 >= np.rad2deg(ICON_grid.clat.min().values) + offset) &
                (s5p_data.latitude
                 <= np.rad2deg(ICON_grid.clat.max().values) - offset),
                drop=True).where(s5p_data.xco2_quality_flag == 0, drop=True)
        except Exception:
            logging.info(
                f"No observations remain after filtering {file} to ICON grid limits"
            )
            _empty_dataset(file[0].name).to_netcdf(output_file)
            continue

        s5p_out = s5p_data[[
            "latitude", "longitude", "date", "xco2", "xco2_quality_flag",
            "xco2_averaging_kernel", "pressure_levels", "pressure_levels",
            "pressure_weight", "co2_profile_apriori", "xco2_apriori",
            "xco2_uncertainty"
        ]]
        s5p_out = s5p_out.rename({
            "levels": "layers",
            "sounding_id": "soundings",
            "xco2": "obs",
            "xco2_quality_flag": "quality_flag",
            "xco2_averaging_kernel": "averaging_kernel",
            "pressure_weight": "pressure_weighting_function",
            "co2_profile_apriori": "prior_profile",
            "xco2_apriori": "prior",
            "xco2_uncertainty": "uncertainty"
        })
        s5p_out["pressure_levels"][:] = s5p_out.pressure_levels[:, ::-1].values
        s5p_out["pressure_weighting_function"][:] = s5p_out.pressure_weighting_function[:, ::-1].values
        s5p_out["prior_profile"][:] = s5p_out.prior_profile[:, ::-1].values
        s5p_out["surface_pressure"] = s5p_out.pressure_levels[:, 0]
        s5p_out.attrs.update({
            'creation_date': str(datetime.now()),
            'author': 'Processing Chain',
            'level_def': 'pressure_boundaries',
            'retrieval_id': file[0].name
        })
        logging.info(
            f"Finished processing and stored obs for date {day.strftime('%Y-%m-%d')}, output file: {output_file} with {s5p_out.soundings.size} soundings"
        )
        s5p_out.to_netcdf(output_file)
