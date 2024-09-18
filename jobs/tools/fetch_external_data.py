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


def fetch_era5(date, dir2move):
    """Fetch ERA5 data from ECMWF for initial conditions

    Parameters
    ----------
    date : initial date to fetch

    """

    c = cdsapi.Client()

    # -- CRWC : Specific rain water content              - 75
    # -- CSWC : Specific snow water content              - 76
    # -- T    : Temperature                             - 130
    # -- U    : U component of wind                     - 131
    # -- V    : V component of wind                     - 132
    # -- Q    : Specific humidity                       - 133
    # -- W    : Vertical velocity                       - 135
    # -- CLWC : Specific cloud liquid water content     - 246
    # -- CIWC : Specific cloud ice water content        - 247

    c.retrieve(
        'reanalysis-era5-complete', {
            'class': 'ea',
            'date': date.strftime('%Y-%m-%d'),
            'time': date.strftime('%H:%M:%S'),
            'expver': '1',
            'levelist':
            '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
            'levtype': 'ml',
            'param': '75/76/130/131/132/133/135/246/247',
            'stream': 'oper',
            'type': 'an',
            'grid': '1.0/1.0',
        }, 'era5_ml.grib')

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

    c.retrieve(
        'reanalysis-era5-single-levels', {
            'product_type': 'reanalysis',
            'param':
            '31/32/33/34/39/40/41/42/43/129/134/139/141/170/172/183/198/235/236/238',
            'date': date.strftime('%Y-%m-%d'),
            'time': date.strftime('%H:%M:%S'),
            'grid': '1.0/1.0',
        }, 'era5_surf.grib')

    shutil.move('era5_ml.grib', os.path.join(dir2move, 'era5_ml.grib'))
    shutil.move('era5_surf.grib', os.path.join(dir2move, 'era5_surf.grib'))


def fetch_era5_nudging(date, dir2move):
    """Fetch ERA5 data from ECMWF for global nudging

    Parameters
    ----------
    date : initial date to fetch

    """

    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-complete', {
            'class': 'ea',
            'date': date.strftime('%Y-%m-%d'),
            'time': date.strftime('%H:%M:%S'),
            'expver': '1',
            'levelist':
            '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
            'levtype': 'ml',
            'param': '75/76/130/131/132/133/135/246/247',
            'stream': 'oper',
            'type': 'an',
            'grid': '1.0/1.0',
        }, 'era5_ml_nudging.grib')

    c.retrieve(
        'reanalysis-era5-single-levels', {
            'product_type': 'reanalysis',
            'param': '129/134',
            'date': date.strftime('%Y-%m-%d'),
            'time': date.strftime('%H:%M:%S'),
            'grid': '1.0/1.0',
        }, 'era5_surf_nudging.grib')

    shutil.move('era5_ml_nudging.grib',
                os.path.join(dir2move, 'era5_ml_nudging.grib'))
    shutil.move('era5_surf_nudging.grib',
                os.path.join(dir2move, 'era5_surf_nudging.grib'))


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

    c = cdsapi.Client()

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


def fetch_OCO2(starttime, endtime, minlon, maxlon, minlat, maxlat, output_folder, product="OCO2_L2_Lite_FP_11r"):

    # hmm. Not currently working. The data is there, https://oco2.gesdisc.eosdis.nasa.gov/data/OCO2_DATA/OCO2_L2_Lite_FP.11.1r/2020/ 
    # but GES DISC doesn't currently have it anymore...
    
    # Set the product (based on the list above!) and other output settings
    product = product # Standard
    begTime = f'{starttime.strftime("%Y-%m-%d")}T00:00:00.000Z'
    endTime = f'{endtime.strftime("%Y-%m-%d")}T23:59:59.999Z'

    # Create a urllib PoolManager instance to make requests.
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',ca_certs=certifi.where())

    # Set the URL for the GES DISC subset service endpoint
    svcurl = 'https://disc.gsfc.nasa.gov/service/subset/jsonwsp'

    # This method POSTs formatted JSON WSP requests to the GES DISC endpoint URL
    # It is created for convenience since this task will be repeated more than once
    def get_http_data(request):
        hdrs = {'Content-Type': 'application/json',
                'Accept'      : 'application/json'}
        data = json.dumps(request)       
        r = http.request('POST', svcurl, body=data, headers=hdrs)
        response = json.loads(r.data)   
        # Check for errors
        if response['type'] == 'jsonwsp/fault' :
            print('API Error: faulty request')
            sys.exit(1)
        return response

    # Construct JSON WSP request for API method: subset
    subset_request = {
        'methodname': 'subset',
        'type': 'jsonwsp/request',
        'version': '1.0',
        'args': {
            'role'  : 'subset',
            'start' : begTime,
            'end'   : endTime,
            'box'   : [minlon, minlat, maxlon, maxlat],  
            'crop'  : False,  
            'data'  : [{'datasetId': product}]
        }
    }

    print("still here...1")

    # Submit the subset request to the GES DISC Server
    response = get_http_data(subset_request)

    print("and even here?")

    # Report the JobID and initial status
    myJobId = response['result']['jobId']
    print('Job ID: '+myJobId)
    print('Job status: '+response['result']['Status'])

    # Construct JSON WSP request for API method: GetStatus
    status_request = {
        'methodname': 'GetStatus',
        'version': '1.0',
        'type': 'jsonwsp/request',
        'args': {'jobId': myJobId}
    }

    # Check on the job status after a brief nap
    while response['result']['Status'] in ['Accepted', 'Running']:
        sleep(5)
        response = get_http_data(status_request)
        status  = response['result']['Status']
        percent = response['result']['PercentCompleted']
        print ('Job status: %s (%d%c complete)' % (status,percent,'%'))

    if response['result']['Status'] == 'Succeeded' :
        print ('Job Finished:  %s' % response['result']['message'])
    else : 
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
    while count < total :
        results_request['args']['startIndex'] += batchsize 
        response = get_http_data(results_request) 
        count = count + response['result']['itemsPerPage']
        results.extend(response['result']['items'])
            
    # Check on the bookkeeping
    print('Retrieved %d out of %d expected items' % (len(results), total))

    # Sort the results into documents and URLs
    docs = []
    urls = []
    for item in results :
        try:
            if item['start'] and item['end'] : urls.append(item) 
        except:
            docs.append(item)

    # Print out the documentation links, but do not download them
    print('\nDocumentation:')
    for item in docs : print(item['label']+': '+item['link'])
        
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
            f = open(outfn,'wb')
            f.write(result.content)
            f.close()
            print(outfn, URL)
        except:
            print('Error! Status code is %d for this URL:\n%s' % (result.status.code,URL))
            print('Help for downloading data is at https://disc.gsfc.nasa.gov/data-access')
            
    print('Finished')
