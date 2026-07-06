import numpy as np
import xarray as xr
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import dask


MONTH_ABBR = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
MONTH_NUM  = {abbr: i + 1 for i, abbr in enumerate(MONTH_ABBR)}


def _parse_season(season_str):
    """Return (month_nums, anchor_month, crosses_year) for a season string."""
    season_str = season_str.upper()
    
    # Direct mapping for standard 3-letter season abbreviations
    SEASON_MAP = {
        'DJF': [12, 1, 2], 'JFM': [1, 2, 3], 'FMA': [2, 3, 4],
        'MAM': [3, 4, 5], 'AMJ': [4, 5, 6], 'MJJ': [5, 6, 7],
        'JJA': [6, 7, 8], 'JAS': [7, 8, 9], 'ASO': [8, 9, 10],
        'SON': [9, 10, 11], 'OND': [10, 11, 12], 'NDJ': [11, 12, 1]
    }
    
    if season_str in SEASON_MAP:
        months = SEASON_MAP[season_str]
    else:
        # Fallback just in case a string like 'DECJANFEB' is passed
        months = [MONTH_NUM[season_str[i*3:(i+1)*3]] for i in range(len(season_str) // 3)]
        
    crosses_year = any(months[i] > months[i + 1] for i in range(len(months) - 1))
    anchor_month = months[-1]
    
    return months, anchor_month, crosses_year


def _assign_season_year(time_da, anchor_month, crosses_year):
    """
    Assign each timestamp its season-year label.
    Season year = calendar year of the last (anchor) month in the season.
    For cross-year seasons, months whose number exceeds the anchor belong
    to the season whose anchor falls in the next calendar year → add 1.
    """
    cal_year  = time_da['time.year']
    cal_month = time_da['time.month']

    if crosses_year:
        shift = (cal_month > anchor_month).astype(int)
        return (cal_year + shift).rename('year')
    else:
        return cal_year.rename('year')


# ---------------------------------------------------------------------------
# format_data  (previously hardcoded to DJF)
# ---------------------------------------------------------------------------

def format_data(data, season='DJF'):
    """
    Select `season` months, assign the correct season-year label, and return
    a days-in-month-weighted mean DataArray indexed by 'year'.

    Convention: the season year equals the calendar year of the *last* month
    in the season string.
    """
    
    with dask.config.set(**{'array.slicing.split_large_chunks': False}):
    
        if 'lon' not in data.coords:
            data = data.rename({'longitude': 'lon', 'latitude': 'lat'})
    
        if data['lon'].max() > 180:
            data = data.assign_coords(lon=((data.lon + 180) % 360) - 180)
            data = data.sortby('lon')
        
        month_nums, anchor_month, crosses_year = _parse_season(season)
    
        data_season = data.sel(time=data['time'].dt.month.isin(month_nums))
    
        year_coord = _assign_season_year(data_season, anchor_month, crosses_year)
        data_season = data_season.assign_coords(year=year_coord)
    
        days_in_month = data_season['time'].dt.days_in_month.assign_coords(year=year_coord)

        weighted_sum = (data_season * days_in_month).groupby('year')

        weighted_sum = weighted_sum.sum(dim='time')

        total_weights = days_in_month.groupby('year').sum(dim='time')

    
    return weighted_sum / total_weights


# ---------------------------------------------------------------------------
# Correlation and significance testing (unchanged)
# ---------------------------------------------------------------------------

def corr_significance(met, enso):
    common_years = np.intersect1d(met['year'], enso.dropna(dim='year')['year'])
    met = met.sel(year=common_years)
    enso = enso.sel(year=common_years)
    n = len(common_years)

    # vectorized, no per-pixel Python calls
    met_anom = met - met.mean('year')
    enso_anom = enso # - enso.mean('year') # enso index already are anomalies
    print(enso.mean('year'))

    cov = (met_anom * enso_anom).sum('year')
    corr = cov / np.sqrt((met_anom**2).sum('year') * (enso_anom**2).sum('year'))

    t_stat = corr * np.sqrt((n - 2) / (1 - corr**2))

    pval = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=n - 2))

    pval = xr.DataArray(pval, coords=corr.coords, dims=corr.dims)

    return corr, pval, n


def BH_test(pvals, alpha=0.05):
    """Benjamini-Hochberg FDR correction applied over the full spatial grid."""
    pval_data = pvals.values.flatten()

    pval_nan  = np.where(np.isnan(pval_data))
    pval_mask = pval_data[~np.isnan(pval_data)]

    _, adjusted_pval, _, _ = multipletests(pval_mask, alpha=alpha, method='fdr_bh')

    fix_pval = np.empty_like(pval_data)
    fix_pval[pval_nan]              = np.nan
    fix_pval[~np.isnan(pval_data)]  = adjusted_pval

    BH_pval = fix_pval.reshape(pvals.shape)
    BH_pval = xr.DataArray(BH_pval, dims=pvals.dims, coords=pvals.coords)

    return BH_pval
