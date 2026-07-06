import pandas as pd
import xarray as xr
import os
import numpy as np
from scipy.signal import find_peaks


from supp_functions.calc_functions import format_data

# ---------------------------------------------------------------------------
# Wavelength identification (unchanged)
# ---------------------------------------------------------------------------

def identify_wavelengths(ds):
    #ds_subset = ds.sel(latitude=slice(31, 36), longitude=slice(267, 278)).load()
    ds_mean = ds.sel(latitude=33, longitude=270, method='nearest')
    #ds_mean = ds_subset.mean(dim=['latitude', 'longitude'])
    #print('data averaged')
    print('data subsetted')

    time_periods = {
        'Original': [],
        'Noise': [],
        'Synoptic': [],
        'Subseasonal': [],
        'Interannual': [],
        'Decadal': [],
        'Secular': []
    }

    time_periods['Original'].append(0)
    time_periods['Noise'].append(1)
    start_mode = 2
    factor = 30

    # check if daily
    if np.median(np.diff(ds_mean.time.values).astype('timedelta64[D]')) == 1:
        time_periods['Synoptic'].append(2)
        start_mode += 1
        factor = 1

    for mode in range(start_mode, ds_mean.shape[0]):
        print(mode)
        arr = ds_mean[mode, :]
        arr_series = pd.Series(arr)

        prominence = np.max(arr_series) * 0.005
        peaks, _ = find_peaks(arr_series, prominence=prominence)

        if np.all(np.diff(arr_series) >= 0) or np.all(np.diff(arr_series) <= 0):
            time_periods['Secular'].append(mode)

        elif len(peaks) <= 1:
            time_periods['Decadal'].append(mode)

        else:
            diffs = np.diff(peaks)
            min_wavelength = np.percentile(diffs, 25) * factor
            max_wavelength = np.percentile(diffs, 75) * factor

            if min_wavelength > 10 and max_wavelength <= 90:
                time_periods['Subseasonal'].append(mode)
            elif min_wavelength > 90 and max_wavelength < 3650:
                time_periods['Interannual'].append(mode)
            elif min_wavelength >= 3650:
                time_periods['Decadal'].append(mode)
            else:
                print(f'Mode {mode} does not fit into any group')

    return time_periods


# ---------------------------------------------------------------------------
# ENSO data reader
# ---------------------------------------------------------------------------


def read_enso_data(sheet_name, months='DJF', interannual=False):
    """
    Read the ENSO index column for `months` season from the Excel file.
    The 'year' column in the Excel file must use the same anchor-month
    convention as process_MET_VAR (i.e. year of the last month in the season).
    """
    
    if interannual:
        print('using interannual indices...')
        enso_file = "/unity/f1/cemerson/data/ENSO_filter/ENSO_Interannual.csv"
        enso = pd.read_csv(enso_file)
        
        enso['time'] = pd.to_datetime(enso['time'])
        enso = enso.set_index('time')
        
        enso = enso[sheet_name]
        time_name = 'time'
        
    
    else:
        enso_file = '/unity/f1/cemerson/data/ENSO_indexes.xlsx'
        enso_df = pd.read_excel(enso_file, sheet_name=sheet_name)
        enso = enso_df[months]
        
        enso.index = enso_df['year']
        time_name = 'year'


    enso_da = xr.DataArray(
        enso,
        coords=[enso.index],
        dims=[time_name],
        name=sheet_name
    )

    mask = (enso_da != 999)   # filter JMA 999 fill value
    if interannual:
        enso_filtered = enso_da.sel(time=mask)
    else:
        enso_filtered = enso_da.sel(year=mask)
    
    if interannual:
        ds_djf = enso_da.sel(time=enso_da['time'].dt.month.isin([12, 1, 2]))
    
        shifted_year = (ds_djf['time.year'] + (ds_djf['time.month'] == 12).astype(int)).rename('year')
        ds_djf.coords['year'] = shifted_year

        weighted_djf = ds_djf.groupby('year').apply(
          lambda group: group.weighted(group.time.dt.days_in_month).mean(dim='time'))
          
        return weighted_djf
        
    else:
        return enso_filtered
