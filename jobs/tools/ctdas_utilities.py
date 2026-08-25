import xarray as xr
import numpy as np


def create_prior_all_ones(output_path,
                          nensembles,
                          ncats,
                          nregs,
                          propagate_bg=False):
    """
    Create a dataset of initial lambdas (all ones) for testing.
    """
    nensembles = nensembles + 1 if propagate_bg else nensembles
    arr = np.ones((nensembles, nregs, ncats, 1), dtype=np.float32)
    arr[-1, :, :, :] = 0 if propagate_bg else 1
    data = xr.DataArray(arr, dims=['ens', 'reg', 'cat', 'tracer'])
    ds = xr.Dataset({'lambda': data})
    try:
        ds.to_netcdf(output_path)
    except:
        print("File currently open. Please close the file and try again.")
    print(f"Prior all ones saved to {output_path}")


def create_prior_all_zeros(output_path, nensembles, ncats, nregs):
    """
    Create a dataset of initial lambdas (all zeros) for testing.
    """
    # Create an array of zeros for all ensembles/regions/categories
    arr = np.zeros((nensembles, nregs, ncats, 1), dtype=np.float32)
    # Add one extra member that is all ones for all regions/categories
    arr = np.vstack((arr, np.ones((1, nregs, ncats, 1), dtype=np.float32)))
    data = xr.DataArray(arr, dims=['ens', 'reg', 'cat', 'tracer'])
    ds = xr.Dataset({'lambda': data})
    try:
        ds.to_netcdf(output_path)
    except:
        print("File currently open. Please close the file and try again.")
    print(f"Prior all zeros saved to {output_path}")


def create_boundary_prior_all_ones(output_path,
                                   n_bg_ens,
                                   nensembles,
                                   propagate_bg=False,
                                   author='Processing Chain',
                                   email=None):
    """
    Create boundary lambdas dataset and save to NetCDF.
    """
    nensembles = nensembles + 1 if propagate_bg else nensembles
    lambdas = np.ones((nensembles, n_bg_ens), dtype=np.float32)
    attrs = {'author': author}
    if email:
        attrs['email'] = email
    ds_lambdas = xr.Dataset(data_vars={'lambda': (['ens', 'reg'], lambdas)},
                            coords={
                                'ens': (['ens'], np.arange(nensembles)),
                                'reg': (['reg'], np.arange(n_bg_ens))
                            },
                            attrs=attrs)
    try:
        ds_lambdas.to_netcdf(output_path)
    except:
        print("File currently open. Please close the file and try again.")
    print(f"Boundary lambdas saved to {output_path}")


def create_boundary_prior_separate(output_path,
                                   n_bg_ens,
                                   author='Processing Chain',
                                   email=None):
    """
    Create boundary lambdas dataset and save to NetCDF.
    """
    # One BG region per ensemble member
    lambdas = np.identity(n_bg_ens, dtype=np.float32)
    # Add one extra member that is all ones for the ensemble member
    lambdas = np.vstack((lambdas, np.ones((1, n_bg_ens), dtype=np.float32)))
    attrs = {'author': author}
    if email:
        attrs['email'] = email
    ds_lambdas = xr.Dataset(data_vars={'lambda': (['ens', 'reg'], lambdas)},
                            coords={
                                'ens': (['ens'], np.arange(n_bg_ens + 1)),
                                'reg': (['reg'], np.arange(n_bg_ens))
                            },
                            attrs=attrs)
    try:
        ds_lambdas.to_netcdf(output_path)
    except:
        print("File currently open. Please close the file and try again.")
    print(f"Boundary-separated lambdas saved to {output_path}")
