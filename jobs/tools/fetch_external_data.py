import os
import shutil
import cdsapi
import zipfile
import logging
import xarray as xr
from icoscp.dobj import Dobj
from icoscp.sparql.runsparql import RunSparql
from icoscp_core.icos import bootstrap
from icoscp import cpauth
import numpy as np
import sys
import json
import datetime
import certifi
import urllib3
import requests
from time import sleep
from subprocess import Popen
from datetime import datetime, timedelta


def fetch_CDS(product, date, levels, params, resolution, area, outloc):
    # Obtain CDS authentification from file
    url_cmd = f"grep 'cds' ~/.cdsapirc"
    url = os.popen(url_cmd).read().strip().split(": ")[1]
    key_cmd = f"sed -n '/cds/ {{n;p}}' ~/.cdsapirc"
    key = os.popen(key_cmd).read().strip().split(": ")[1]
    c = cdsapi.Client(url=url, key=key)

    # Set temporal choices. ERA5 data on disk uses lists [2018-01-01, 2018-01-02, etc] while ERA5-complete uses strings with / as the separator
    if isinstance(date, datetime):
        datestr = date.strftime('%Y-%m-%d')
        timestr = date.strftime('%H:%M')
    elif isinstance(date, list):
        datestr = sorted({dt.date().strftime("%Y-%m-%d") for dt in date})
        datestr = datestr if product == 'reanalysis-era5-single-levels' else '/'.join(
            map(str, datestr))
        timestr = sorted({dt.time().strftime("%H:%M") for dt in date})
        timestr = timestr if product == 'reanalysis-era5-single-levels' else '/'.join(
            map(str, timestr))
    else:
        raise TypeError(
            f"Expected a datetime or list, but got {type(date).__name__}.")

    # Set level choices
    if isinstance(levels, str):
        levelstr = levels
    elif isinstance(levels, list):
        levelstr = '/'.join(map(str, levels))
    elif levels is None:
        pass
    else:
        raise TypeError(
            f"Expected a string or list, but got {type(levels).__name__}.")

    # Set parameters
    if isinstance(params, str):
        paramstr = params
    elif isinstance(params, list):
        paramstr = '/'.join(map(str, params))
    else:
        raise TypeError(
            f"Expected a string or list, but got {type(params).__name__}.")

    c.retrieve(
        product, {
            'date':
            datestr,
            'time':
            timestr,
            'param':
            paramstr,
            'grid':
            f'{resolution}/{resolution}',
            **({'area' : area} if area is not None else {}),
            **({
                'class': 'ea',
                'type': 'an',
                'stream': 'oper',
                'levelist': levelstr,
                'levtype': 'ml',
                'expver': '1'
            } if product == 'reanalysis-era5-complete' else {}),
            **({
                'product_type': 'reanalysis'
            } if product == 'reanalysis-era5-single-levels' else {}),
        }, outloc)


def fetch_era5(date, dir2move, resolution=1.0, area=None):
    if isinstance(date, list):
        outfile_3D = dir2move / f"era5_ml_{date[0].strftime('%Y-%m-%d')}_{date[-1].strftime('%Y-%m-%d')}.grib"
    else:
        outfile_3D = dir2move / f"era5_ml_{date.strftime('%Y-%m-%d')}.grib"
    if not os.path.isfile(outfile_3D):
        # -- CRWC : Specific rain water content              - 75
        # -- CSWC : Specific snow water content              - 76
        # -- T    : Temperature                             - 130
        # -- U    : U component of wind                     - 131
        # -- V    : V component of wind                     - 132
        # -- Q    : Specific humidity                       - 133
        # -- W    : Vertical velocity                       - 135
        # -- CLWC : Specific cloud liquid water content     - 246
        # -- CIWC : Specific cloud ice water content        - 247
        fetch_CDS('reanalysis-era5-complete', date, '1/to/137',
                  [75, 76, 130, 131, 132, 133, 135, 246, 247], resolution, area,
                  outfile_3D)

    if isinstance(date, list):
        outfile_surface = dir2move / f"era5_surf_{date[0].strftime('%Y-%m-%d')}_{date[-1].strftime('%Y-%m-%d')}.grib"
    else:
        outfile_surface = dir2move / f"era5_surf_{date.strftime('%Y-%m-%d')}.grib"
    if not os.path.isfile(outfile_surface):
        # -- CI   : Sea Ice Cover                   - 31
        # -- ASN  : Snow albedo                     - 32
        # -- RSN  : Snow density                    - 33
        # -- SST  : Sea Surface Temperature         - 34
        # -- SWV1 : Volumetric soil water layer 1   - 39
        # -- SWV2 : Volumetric soil water layer 2   - 40
        # -- SWV3 : Volumetric soil water layer 3   - 41
        # -- SWV4 : Volumetric soil water layer 4   - 42
        # -- SLT  : Soil type                       - 43
        # -- Z    : Geopotential                   - 129
        # -- SP   : Surface pressure               - 134
        # -- STL1 : Soil temperature level 1       - 139
        # -- SD   : Snow depth                     - 141
        # -- STL2 : Soil temperature level 2       - 170
        # -- LSM  : Land-Sea Mask                  - 172
        # -- STL3 : Soil temperature level 3       - 183
        # -- SRC  : Skin reservoir content         - 198
        # -- SKT  : Skin Temperature               - 235
        # -- STL4 : Soil temperature level 4       - 236
        # -- TSN  : Temperature of snow layer      - 238
        fetch_CDS('reanalysis-era5-single-levels', date, None, [
            31, 32, 33, 34, 39, 40, 41, 42, 43, 129, 134, 139, 141, 170, 172,
            183, 198, 235, 236, 238
        ], resolution, area, outfile_surface)
    
    return outfile_3D, outfile_surface


def fetch_era5_nudging(date, dir2move, resolution=1.0, area=None):
    """Fetch ERA5 data from ECMWF for global nudging

    Parameters
    ----------
    date : initial date to fetch

    """
    if isinstance(date, list):
        outfile_3D = dir2move / f"era5_ml_nudging_{date[0].strftime('%Y-%m-%d')}_{date[-1].strftime('%Y-%m-%d')}.grib"
    else:
        outfile_3D = dir2move / f"era5_ml_nudging_{date.strftime('%Y-%m-%d')}.grib"
    if not os.path.isfile(outfile_3D):
        fetch_CDS('reanalysis-era5-complete', date, '1/to/137',
                  [75, 76, 130, 131, 132, 133, 135, 246, 247], resolution, area,
                  outfile_3D)

    if isinstance(date, list):
        outfile_surface = dir2move / f"era5_surf_nudging_{date[0].strftime('%Y-%m-%d')}_{date[-1].strftime('%Y-%m-%d')}.grib"
    else:
        outfile_surface = dir2move / f"era5_surf_nudging_{date.strftime('%Y-%m-%d')}.grib"
    if not os.path.isfile(outfile_surface):
        fetch_CDS('reanalysis-era5-single-levels', date, None, [
            129, 134
        ], resolution, area, outfile_surface)

    return outfile_3D, outfile_surface


def fetch_CAMS_CO2(date, dir2move):
    """Fetch CAMS CO2 data from ECMWF for initial and boundary conditions

    Parameters
    ----------
    date : initial date to fetch a year's worth of data

    """

    # Set a temporary destionation
    tmpdir = os.path.join(os.getenv('SCRATCH'), 'CAMS_i')
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    url_cmd = f"grep 'ads' ~/.cdsapirc"
    url = os.popen(url_cmd).read().strip().split(": ")[1]
    key_cmd = f"sed -n '/ads/ {{n;p}}' ~/.cdsapirc"
    key = os.popen(key_cmd).read().strip().split(": ")[1]
    c = cdsapi.Client(url=url, key=key)

    download = os.path.join(tmpdir, f'cams_GHG_{date.strftime("%Y")}.zip')
    if not os.path.isfile(download):
        c.retrieve(
            'cams-global-greenhouse-gas-inversion', {
                'variable':
                'carbon_dioxide',
                'quantity':
                'concentration',
                'input_observations':
                'surface',
                'time_aggregation':
                'instantaneous',
                'version':
                'latest',
                'year':
                date.strftime('%Y'),
                'month': [
                    '01',
                    '02',
                    '03',
                    '04',
                    '05',
                    '06',
                    '07',
                    '08',
                    '09',
                    '10',
                    '11',
                    '12',
                ],
                'format':
                'zip',
            }, download)
        logging.info(f'downloaded the CAMS data!')
    else:
        logging.info(f'File already downloaded and present at {download}')

    # --- Extract the zip file
    with zipfile.ZipFile(download) as zf:
        for member in zf.infolist():
            if not os.path.isfile(os.path.join(tmpdir, member.filename)):
                try:
                    zf.extract(member, tmpdir)
                except zipfile.error as e:
                    pass

    # --- Output files to folder
    with zipfile.ZipFile(download) as zf:
        for member in zf.infolist():
            filename = os.path.join(tmpdir, member.filename)
            logging.info("Writing out CAMS data to file")
            ds_CAMS = xr.open_dataset(filename)
            for time in ds_CAMS.time:
                outpath = os.path.join(
                    dir2move, 'cams_egg4_' +
                    ds_CAMS.sel(time=time).time.dt.strftime('%Y%m%d%H').values
                    + '.nc')
                if not os.path.isfile(outpath):
                    ds_out = ds_CAMS.where(ds_CAMS.time == time,
                                           drop=True).squeeze()
                    ds_out.to_netcdf(outpath)


def fetch_ICOS_data(cookie_token,
                    query_type='any',
                    start_date='01-01-2022',
                    end_date='31-12-2022',
                    save_path='',
                    species=['co', 'co2', 'ch4']):
    '''
    This script starts a SPARQL query for downloading ICOS-CP data. The query is based on searching at the ICOS-CP
    (e.g., https://data.icos-cp.eu/portal/#%7B%22filterCategories%22%3A%7B%22variable%22%3A%5B%22http%3A%2F%2Fmeta.icos-cp.eu%2Fresources%2Fcpmeta%2Fco2atcMoleFrac%22%5D%7D%2C%22filterTemporal%22%3A%7B%22df%22%3A%222017-12-31%22%2C%22dt%22%3A%222018-12-30%22%7D%7D)
    and then clicking the well-hidden SPARQL query button (situated right of "Data objects 1 to 20 of 167", consisting of an arrow.)

    cookie_token    str    cpauthToken=WzE3M....
    query_type      str    [release, growing, any] correspond to the different file products at the ICOS-CP
    start_date      str    dd-mm-yyyy
    end_date        str    dd-mm-yyyy
    save_path       str    e.g., /scratch/snx/[user]/ICOS_data/year/
    species         list   can be ['co', 'co2', 'ch4'] or any subset thereof
    '''
    meta, data = bootstrap.fromCookieToken(cookie_token)
    cpauth.init_by(data.auth)
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

    for d in result.data()['dobj']:
        obj = Dobj(d).data

        shape = np.shape(obj)

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
        try:
            ds = ds.sel(index=slice(
                datetime.strptime(start_date, '%d-%m-%Y').date().strftime(
                    '%Y-%m-%d'),
                datetime.strptime(end_date, '%d-%m-%Y').date().strftime(
                    '%Y-%m-%d')))
        except:
            print('failure!')
            print(ds.index)
            print(f"Not doing {Dobj(d).station['id']}, then...?")
            break
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
        name = 'ICOS_obs_' + str(specie)[2:-2] + '_' + query_type + '_' + str(
            Dobj(d).station['id']) + '_' + str(
                Dobj(d).meta['specificInfo']['acquisition']
                ['samplingHeight']) + '_' + start_date + '_' + end_date + '.nc'
        ds.to_netcdf(os.path.join(save_path, name))


def fetch_OCO2(starttime,
               endtime,
               minlon,
               maxlon,
               minlat,
               maxlat,
               output_folder,
               product="OCO2_L2_Lite_FP_11r"):

    # hmm. Not currently working. The data is there, https://oco2.gesdisc.eosdis.nasa.gov/data/OCO2_DATA/OCO2_L2_Lite_FP.11.1r/2020/
    # but GES DISC doesn't currently have it anymore...

    # Set the product (based on the list above!) and other output settings
    product = product  # Standard
    begTime = f'{starttime.strftime("%Y-%m-%d")}T00:00:00.000Z'
    endTime = f'{endtime.strftime("%Y-%m-%d")}T23:59:59.999Z'

    # Create a urllib PoolManager instance to make requests.
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',
                               ca_certs=certifi.where())

    # Set the URL for the GES DISC subset service endpoint
    svcurl = 'https://disc.gsfc.nasa.gov/service/subset/jsonwsp'

    # This method POSTs formatted JSON WSP requests to the GES DISC endpoint URL
    # It is created for convenience since this task will be repeated more than once
    def get_http_data(request):
        hdrs = {
            'Content-Type': 'application/json',
            'Accept': 'application/json'
        }
        data = json.dumps(request)
        r = http.request('POST', svcurl, body=data, headers=hdrs)
        response = json.loads(r.data)
        # Check for errors
        if response['type'] == 'jsonwsp/fault':
            print('API Error: faulty request')
            sys.exit(1)
        return response

    # Construct JSON WSP request for API method: subset
    subset_request = {
        'methodname': 'subset',
        'type': 'jsonwsp/request',
        'version': '1.0',
        'args': {
            'role': 'subset',
            'start': begTime,
            'end': endTime,
            'box': [minlon, minlat, maxlon, maxlat],
            'crop': False,
            'data': [{
                'datasetId': product
            }]
        }
    }

    print("still here...1")

    # Submit the subset request to the GES DISC Server
    response = get_http_data(subset_request)

    print("and even here?")

    # Report the JobID and initial status
    myJobId = response['result']['jobId']
    print('Job ID: ' + myJobId)
    print('Job status: ' + response['result']['Status'])

    # Construct JSON WSP request for API method: GetStatus
    status_request = {
        'methodname': 'GetStatus',
        'version': '1.0',
        'type': 'jsonwsp/request',
        'args': {
            'jobId': myJobId
        }
    }

    # Check on the job status after a brief nap
    while response['result']['Status'] in ['Accepted', 'Running']:
        sleep(5)
        response = get_http_data(status_request)
        status = response['result']['Status']
        percent = response['result']['PercentCompleted']
        print('Job status: %s (%d%c complete)' % (status, percent, '%'))

    if response['result']['Status'] == 'Succeeded':
        print('Job Finished:  %s' % response['result']['message'])
    else:
        print('Job Failed: %s' % response['fault']['code'])
        sys.exit(1)

    # Construct JSON WSP request for API method: GetResult
    batchsize = 20
    results_request = {
        'methodname': 'GetResult',
        'version': '1.0',
        'type': 'jsonwsp/request',
        'args': {
            'jobId': myJobId,
            'count': batchsize,
            'startIndex': 0
        }
    }

    # Retrieve the results in JSON in multiple batches
    # Initialize variables, then submit the first GetResults request
    # Add the results from this batch to the list and increment the count
    results = []
    count = 0
    response = get_http_data(results_request)
    count = count + response['result']['itemsPerPage']
    results.extend(response['result']['items'])

    # Increment the startIndex and keep asking for more results until we have them all
    total = response['result']['totalResults']
    while count < total:
        results_request['args']['startIndex'] += batchsize
        response = get_http_data(results_request)
        count = count + response['result']['itemsPerPage']
        results.extend(response['result']['items'])

    # Check on the bookkeeping
    print('Retrieved %d out of %d expected items' % (len(results), total))

    # Sort the results into documents and URLs
    docs = []
    urls = []
    for item in results:
        try:
            if item['start'] and item['end']: urls.append(item)
        except:
            docs.append(item)

    # Print out the documentation links, but do not download them
    print('\nDocumentation:')
    for item in docs:
        print(item['label'] + ': ' + item['link'])

    # Use the requests library to submit the HTTP_Services URLs and write out the results.
    print('\nHTTP_services output:')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    for item in urls:
        outfn = output_folder + '/' + item['label']
        if os.path.isfile(outfn):
            continue

        URL = item['link']
        result = requests.get(URL)
        try:
            result.raise_for_status()
            f = open(outfn, 'wb')
            f.write(result.content)
            f.close()
            print(outfn, URL)
        except:
            print('Error! Status code is %d for this URL:\n%s' %
                  (result.status.code, URL))
            print(
                'Help for downloading data is at https://disc.gsfc.nasa.gov/data-access'
            )

    print('Finished')

def process_OCO2():
    ######### Some messages ######### 
    print('=============================================================================') 
    print(' Pre-processing Observation product, readable for CTDAS-ICON.' ) 
    print(' Data will be filtered base on a given ICON domain.' ) 
    print(' David Ho, MPI-BGC Jena' ) 
    print('=============================================================================') 
    print('') 
    print('Loading neccessary packages...') 
    ## Import 
    import numpy as np 
    import pandas as pd 
    import xarray as xr 
    import glob 
    from netCDF4 import Dataset 
    import datetime 
    import time as TIME 
    import warnings 
    import os
    warnings.filterwarnings("ignore") 
    print('') 
    #-- retrieve start time 
    t1 = TIME.time() 
    ######### Output path ########### 
    nc_out = '//scratch/snx3000/ekoene/OCO-2_filtered/' 
    if not os.path.exists(nc_out):
        os.makedirs(nc_out)
        print(f"Output folder '{nc_out}' created successfully.")
    else:
        print(f"Output folder '{nc_out}' already exists.")
    ######### Time control ########### 
    Year = 2018 
    for month in range(1,13):
        if month in [4, 6, 9, 11]: 
            daymax = 30 
        elif month == 2: 
            daymax = 28 
        else: 
            daymax = 31 

        ndays = np.arange(1, daymax+1) # 1st~31th 
        ######### Observation ########### 
        file_list = sorted( glob.glob('/scratch/snx3000/ekoene/OCO-2/OCO2_L2_Lite_FP.11r:oco2_LtCO2_*') ) 
        if len(file_list) == 0: 
            raise ValueError("File list is empty, stopping here!") 

        ########## ICON grid ############ 
        mainpath = '/users/ekoene/CTDAS_inputs/' 
        grid_file = mainpath + '/icon_europe_DOM01.nc' 
        ICON_GRID = xr.open_dataset(grid_file) 
        # Convert an array of size 1 to its scalar equivalent. 
        lon_min = np.min(ICON_GRID.clon.values) 
        lon_max = np.max(ICON_GRID.clon.values) 
        lat_min = np.min(ICON_GRID.clat.values) 
        lat_max = np.max(ICON_GRID.clat.values) 
        print('ICON grid extends:') 
        print('Longitude min. %7.4f, max. %7.4f' % (np.rad2deg(lon_min),np.rad2deg(lon_max)) ) 
        print('Latitude min. %7.4f, max. %7.4f' % (np.rad2deg(lat_min),np.rad2deg(lat_max)) ) 
        print('') 
        ########## Set bounds to filter ########## 
        offset = 1.2 
        sub_lon_min = np.rad2deg(lon_min) + offset 
        sub_lon_max = np.rad2deg(lon_max) - offset 
        sub_lat_min = np.rad2deg(lat_min) + offset 
        sub_lat_max = np.rad2deg(lat_max) - offset 
        print('To avoid cells at the domain boundary, subtracting: %s degree.' %offset) 
        print('Filtered extends:') 
        print('Longitude min. %7.4f, max. %7.4f' %(sub_lon_min, sub_lon_max) ) 
        print('Latitude min. %7.4f, max. %7.4f' %(sub_lat_min, sub_lat_max) ) 
        print('') 

        ######## Begin Production ############# 
        Total_nobs_before = np.array([]) 
        Total_nobs_after = np.array([]) 
        for day in ndays: 
            print('Processing: (%s/%s)' %(day, len(ndays)) ) 
            ######### Read data #########  
            try: 
                # Find a file in the file list
                for file_name in file_list:
                    if f"OCO2_L2_Lite_FP.11r:oco2_LtCO2_{str(Year)[2:]}{month:02d}{day:02d}" in file_name:
                        s5p_file = file_name 
                        print('Opening file: %s' %s5p_file) 
                        s5p_data = Dataset(s5p_file) 
            except: 
                print('file %s not found.' %s5p_file) 
                print('Skipping...') 
                print('') 
                continue # Continue to next iteration. 

        ######## Filter base of ICON domain ######## 
            date_list = [] 
            for timestamp in s5p_data['time'][:]: 
                value = datetime.datetime.fromtimestamp(timestamp) 
                date_list.append(value) 
            
            dictionary = { 
            'date_time' : date_list[:], 
            'raw_time' : s5p_data['time'][:], 
            'xco2': s5p_data['xco2'][:], 
            'lat': s5p_data['latitude'][:], 
            'lon': s5p_data['longitude'][:], 
            'qf': s5p_data['xco2_quality_flag'][:], # quality flag 0 = good; 1 = bad.
            } 
            df_pixels = pd.DataFrame(data=dictionary) 

            ## Filter base on ICON domain ## 
            inside_domain_flag = ( ( df_pixels['lon'] > sub_lon_min ) & ( df_pixels['lon'] < sub_lon_max )  & 
                                ( df_pixels['lat'] > sub_lat_min ) & ( df_pixels['lat'] < sub_lat_max )  *
                                ( df_pixels['qf'] == 0) )

            # -- Old hard coded settings: 
            # inside_domain_flag = ( ( df_pixels['lon'] > -20 ) & ( df_pixels['lon'] < 58 ) \ 
            # & ( df_pixels['lat'] > 32 ) & ( df_pixels['lat'] < 69 ) ) 
            ## Get the indexes from data frame ## 
            indexes = df_pixels[inside_domain_flag].index 

            ## Some messages 
            Before = len(s5p_data.variables['xco2'][:]) 
            print('It had %i data' %Before) 
            Total_nobs_before = np.append(Total_nobs_before, Before) 

            After = len(s5p_data.variables['xco2'][indexes]) 
            print('Now has %i' %After) 
            Total_nobs_after = np.append(Total_nobs_after, After) 
            if After == 0:
                print('skipping')
                continue

            ######### Create/Write netCDF ######### 
            _, tail = os.path.split(s5p_file)
            output_path = os.path.join(nc_out, 'OCO2_%04d%02d%02d_ctdas.nc' %(Year, month, day)) 
            ncfile = Dataset( output_path, mode='w', format='NETCDF4' ) 
            print('Writing %s from %s' %(output_path, s5p_file)) 

            ######### Def. attribute ######### 
            ncfile.level_def = 'pressure_boundaries' 
            ncfile.retrieval_id = tail
            ncfile.creator_name = 'Erik Koene (Empa)' 
            ncfile.date_created = str( datetime.datetime.now() ) 
            ######### Create dimension ######### 
            #ncfile.createDimension( 'soundings', s5p_data.dimensions['sounding_dim'].size ) # Select all 
            ncfile.createDimension( 'soundings', s5p_data['xco2'][indexes].size ) # Select the indexes 
            ncfile.createDimension( 'levels', s5p_data.dimensions['levels'].size ) 
            ncfile.createDimension( 'layers', s5p_data.dimensions['levels'].size ) 
            ncfile.createDimension( 'epoch_dimension', s5p_data.dimensions['epoch_dimension'].size )
            ######### Set variables ######### 
            ### Lat/Lon 
            lat = ncfile.createVariable('latitude', np.float32, ('soundings')) 
            lat.units = 'degrees_north' 
            #lat[:] = s5p_data.variables['latitude'][:] 
            lat[:] = s5p_data.variables['latitude'][indexes] 
            lon = ncfile.createVariable('longitude', np.float32, ('soundings')) 
            lon.units = 'degrees_east' 
            #lon[:] = s5p_data.variables['longitude'][:] 
            lon[:] = s5p_data.variables['longitude'][indexes] 
            ### Time 
            date = ncfile.createVariable('date', np.uint32, ('soundings', 'epoch_dimension')) 
            date.units = 'seconds since 1970-01-01 00:00:00' # 
            date.long_name = 'date_time' 
            # Converting... 
            A = np.array([], np.uint32) 
            for timestamp in s5p_data['time'][indexes]: 
                value = datetime.datetime.fromtimestamp(timestamp) 
                time = np.array([value.year, value.month, value.day, value.hour, value.minute, value.second, value.microsecond], np.uint32) 
                A = np.concatenate([A, time], axis=0) 
            B = A.reshape( int(len(A)/7), 7 ) 
            date[:] = B[:] 
            ##### Obs 
            obs = ncfile.createVariable('obs', np.float32, ('soundings')) 
            obs.units = '1e-6 [ppm]' 
            obs.long_name = 'column-averaged dry air mole fraction of atmospheric co2' 
            obs.comment = 'Retrieved column-averaged dry air mole fraction of atmospheric carbon dioxide (XCO2) in ppm for CTDAS' 
            obs[:] = s5p_data.variables['xco2'][indexes] # Keep ppm units
            ### qa flag 
            qa_f = ncfile.createVariable('quality_flag', np.int8, ('soundings')) 
            qa_f.flag_values= '[0, 1]' 
            qa_f.long_name = 'quality flag for the retrieved column-averaged dry air mole fraction of atmospheric methane' 
            qa_f.comment = '0=good, 1=bad' 
            qa_f[:] = s5p_data.variables['xco2_quality_flag'][indexes] 
            ##### avg kernel 
            avg_kernel = ncfile.createVariable('averaging_kernel', np.float32, ('soundings', 'layers')) 
            avg_kernel.units = '1' 
            avg_kernel.long_name = 'xco2 averaging kernel' 
            avg_kernel.comment = 'Represents the altitude sensitivity of the retrieval as a function of pressure. All values represent layer averages within the corresponding pressure levels. Profiles are ordered from surface to top of atmosphere.' 
            #avg_kernel[:] = s5p_data.variables['xch4_averaging_kernel'][:] 
            avg_kernel[:] = s5p_data.variables['xco2_averaging_kernel'][:][indexes] 
            ### surface_pressure 
            psurf = ncfile.createVariable('surface_pressure', np.float32, ('soundings')) 
            psurf.long_name = 'Surface pressure' 
            psurf.comment = 'Sliced from: OCO2_pressure_levels[:, 0] in Python. Pressure levels defined at the same levels as the averaging kernel and a priori profile layers. Levels were ordered from top of atmosphere to surface.' 
            psurf.unit = 'hPa' 
            #psurf[:] = s5p_data.variables['pressure_levels'][:, 0] 
            psurf[:] = s5p_data.variables['pressure_levels'][indexes, -1] 
            ##### pressure_levels 
            pres_lvls = ncfile.createVariable('pressure_levels', np.float32, ('soundings', 'layers')) 
            pres_lvls.long_name = 'Pressure levels' 
            pres_lvls.comment = 'Sliced from: s5p_pressure_levels[:, 1:], Python. Pressure levels define the boundaries of the averaging kernel and a priori profile layers. Levels were ordered from top of atmosphere to surface.' 
            pres_lvls.unit = 'hPa'  
            pres_lvls[:] = s5p_data.variables['pressure_levels'][:, ::-1][indexes] 
            #### pressure_weighting_function 
            pwf = ncfile.createVariable('pressure_weighting_function', np.float32, ('soundings', 'layers')) 
            pwf.long_name = 'Pressure weighting function' 
            pwf.comment = 'Layer dependent weights needed to apply the averaging kernels.' 
            pwf[:] = s5p_data.variables['pressure_weight'][:,::-1][indexes] 
            ### prior_profile 
            prior_profile = ncfile.createVariable('prior_profile', np.float32, ('soundings', 'layers')) 
            prior_profile.units = '1e-6 [ppm]' 
            prior_profile.long_name = 'a priori dry air mole fraction profile of atmospheric CO2' 
            prior_profile.comment = 'A priori dry-air mole fraction profile of atmospheric CO2 in ppm. All values represent layer averages within the corresponding pressure levels. Profiles are ordered from top of atmosphere to the surface.' 
            prior_profile[:] = s5p_data.variables['co2_profile_apriori'][:,::-1][indexes]
            ### prior 
            prior = ncfile.createVariable('prior', np.float32, ('soundings')) 
            prior.units = '1e-6 [ppm]' 
            prior.long_name = 'Prior' 
            prior.comment = 'The a priori CO2 profile uses the same formulation as used for TCCON GGG2020 retrievals' 
            prior[:] = s5p_data.variables["xco2_apriori"][indexes]
            ### uncertainty 
            unc = ncfile.createVariable('uncertainty', np.float32, ('soundings')) 
            unc.units = '1e-6 [ppm]' 
            unc.long_name = '1-sigma uncertainty of the retrieved column-averaged dry air mole fraction of atmospheric carbon dioxide' 
            unc.comment = '1-sigma uncertainty of the retrieved column-averaged dry air mole fraction of atmospheric carbon dioxide (XCO2) in ppm' 
            unc[:] = s5p_data.variables['xco2_uncertainty'][indexes]

            ### Extras 
            # unique sounding_id 
            sounding_id = ncfile.createVariable('sounding_id', np.int64, ('soundings')) 
            sounding_id.comment ='Some numbers unique per observation' 
            # sounding_id[:] = s5p_data.variables['sounding_id'][indexes]
            #s5p_obs = s5p_data.variables['xch4'][:] 
            s5p_obs = s5p_data.variables['xco2'][indexes] 
            to_add = int('%04d%02d%02d0000000' %(Year, month, day)) # datetime + 9 zeros 
            sounding_id[:] = np.arange(len(s5p_obs)) + to_add 
            print('Added %s to sounding id, unique per observation' %to_add) 

            # nobs 
            #nobs = ncfile.createVariable('nobs', np.unit32, ('soundings')) 
            #nobs.comment ='Number of observations' 
            #nobs[:] = 
            #print(ncfile) 
            ncfile.close() 
            print('Done! Closing netcdf, proceeding to the next file.') 
            print('') 

    t2 = TIME.time() 
    print('') 
    print('All done! Wallclock time: %0.3f seconds' % (t2-t1)) 
    sum_before = np.sum(Total_nobs_before) 
    sum_after = np.sum(Total_nobs_after) 
    print('Summary:') 
    print('Original nobs: %i. Filtered nobs: %i.' %(sum_before, sum_after) ) 
