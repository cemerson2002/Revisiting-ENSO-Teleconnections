# produces interpolation sensitivity plot

from supp_functions.preproc_functions import read_enso_data
from supp_functions.calc_functions import format_data, corr_significance, BH_test
from supp_functions.plotting_functions import apply_multipanel_map_styling

import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import regionmask
import scipy.stats as stats
import dask

from dask import compute

# --- Configuration & Setup ---
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.use('Agg')
land = regionmask.defined_regions.natural_earth_v5_0_0.land_110
xr.set_options(use_flox=False)

file_dict = {
    '200HEIGHT': "/unity/f1/cemerson/data/unfiltered_height/ERA5_height_200_full_fixed.nc",
    '500HEIGHT': "/unity/f1/cemerson/data/500_unfiltered_height/ERA5_height_500.nc",
    'PRECIP': "/unity/f1/cemerson/data/surface_precip/daily_CHIRPSv2_full.nc",
    'TEMP': "/unity/f1/cemerson/data/surface_precip/daily_CRU_temp.nc"
    }
    
target_enso = 'ONI' 

def _pearsonr_ufunc(x, y):

    if np.isnan(x).all() or np.isnan(y).all():
        return np.nan, np.nan
        
    with np.errstate(divide='ignore', invalid='ignore'):
        corr, pval = stats.pearsonr(x, y)
    return corr, pval

def steiger_overlapping(r13, r23, r12, n):
    """
    Tests if r(var1, ENSO) differs significantly from r(var2, ENSO).
    r13: native ~ ENSO
    r23: interpolated ~ ENSO
    r12: native ~ interpolated  (the shared-variable correction)
    n:   original sample size (NOT inflated by interpolation)
    """

    # Fisher Z-transform the two correlations of interest
    z13 = np.arctanh(r13)
    z23 = np.arctanh(r23)

    # Average of the two correlations (used in correction)
    r_bar = (r13 + r23) / 2
    
    r_bar = np.clip(r_bar, -0.9999, 0.9999)

    # Meng et al. correction factor
    f = (1 - r12) / (2 * (1 - r_bar**2))
    h = (1 - f * r_bar**2) / (1 - r_bar**2)

    # Standard error
    se = np.sqrt(1 / (n - 3) * (
        1 / (1 - r12) + 
        h * (r13**2 + r23**2 - 2 * r12 * r13 * r23) / 
        (2 * (1 - r_bar**2)**2)
    ))

    # Test statistic
    z_stat = (z13 - z23) / se

    # Two-tailed p-value
    p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))
    
    return p_value

#    return {
#        'r_native_ENSO':       r13,
#        'r_interp_ENSO':       r23,
#        'r_native_interp':     r12,
#        'z_statistic':         round(z_stat, 4),
#        'p_value':             round(p_value, 4),
#        'significant':         p_value < 0.05
#    }

# --- Color & Formatting Setup ---
colors = ['navy','royalblue','lightblue','white','lightcoral','crimson','darkred']
bounds = [-1.0, -0.5, -0.1, -0.02, 0.02, 0.1, 0.5, 1.0]
cmap = ListedColormap(colors)
norm = BoundaryNorm(bounds, cmap.N)

ENSO_da = read_enso_data(target_enso, months='DJF', interannual=True)

lat_1deg = np.arange(-90, 91, 1)
lon_1deg = np.arange(-180, 180, 1)

fig_proj = ccrs.PlateCarree()
data_proj = ccrs.PlateCarree()

row, col = 2, 2
fig_labels = [chr(i) for i in range(ord('a'), ord('a') + row * col)]

matplotlib.rcParams['axes.linewidth'] = 0.5
fig, axes_3d = plt.subplots(ncols=col, nrows=row, subplot_kw={'projection': fig_proj},
        figsize=(col*6, row*3.5), 
        constrained_layout=True, sharex='col', sharey='row')
        
axes = axes_3d.flatten()

lazy_results = []
metadata = []

for i, (var_name, file_name) in enumerate(file_dict.items()):
    print(var_name)
    print(file_name)
    # --- Process Dataset ---
    dataset_ds = xr.open_dataset(file_name, chunks={'time': -1, 'lat': 181, 'lon': 360})
        
    var_name = list(dataset_ds.data_vars)[0]
    dataset_ds = dataset_ds[var_name]
        
    for l in ["level", "lev"]:
        if l in dataset_ds.coords:
            dataset_ds = dataset_ds.isel({l:0})
    
    processed_ds = format_data(dataset_ds)
    
    
    # --- Land Mask Generation ---
    mask = None
    if var_name in ['TEMP', 'PRECIP']:
        mask = land.mask(processed_ds['lon'], processed_ds['lat'])
    
    # --- Correlation in Native Resolution ---
    print('calculating correlations')
    corr_native, _, n = corr_significance(processed_ds, ENSO_da)
    
    # --- Correlation via Interpolation ---
        
    ds_1deg = processed_ds.interp(lat=lat_1deg, lon=lon_1deg, method='linear')
    
    corr_1deg, _, n = corr_significance(ds_1deg, ENSO_da)
    
    corr_restored = corr_1deg.interp(lat=processed_ds['lat'], lon=processed_ds['lon'], method='linear')
    
    if var_name in ['TEMP', 'PRECIP']:
        corr_restored = corr_restored.where(mask == 0)
        corr_native = corr_native.where(mask == 0)
    
    
    # --- Plotting Differences ---
    corr_diff = corr_restored - corr_native
    
    with np.errstate(divide='ignore', invalid='ignore'):

        print('starting calculations')
        # Apply the function lazily over the 'year' core dimension
        processed_ds = processed_ds.chunk({'year': -1})
        ds_1deg  = ds_1deg.chunk({'year': -1})
        
        corr_both, _ = xr.apply_ufunc(
            _pearsonr_ufunc,
            processed_ds,
            ds_1deg.interp(lat=processed_ds['lat'], lon=processed_ds['lon'], method='linear'),
            input_core_dims=[['year'], ['year']],
            output_core_dims=[[], []],
            vectorize=True,
            dask='parallelized',
            output_dtypes=[float, float]
        )
        
        print(corr_both)
        
    
    print("start Steiger's Test")
    # Steiger's Test for overlapping correlations
    
    
    significance = xr.apply_ufunc(
        steiger_overlapping,
        corr_native,
        corr_restored,
        corr_both,   # this is r12
        input_core_dims=[[], [], []],
        kwargs={'n': n},
        vectorize=True,
        dask='parallelized',
        output_dtypes=[float]
    )
    
    print(significance)
    
    lazy_results.append((corr_diff, significance))
    metadata.append(var_name)
    
    

computed = compute(*[
    item for group in lazy_results for item in group
])
print('done computing')

grouped_results = [computed[i:i+2] for i in range(0, len(computed), 2)]

print('beginning plotting')
# ---------------------------
# Plotting (now NumPy-backed)
# ---------------------------

print(np.shape(grouped_results))

for i, (corr_diff, significance) in enumerate(grouped_results):  
    print(corr_diff)
    print(significance)
    print(np.max(np.abs(corr_diff)))
    print(np.max(corr_diff))
    print(np.min(corr_diff))
    
    print(np.max(significance))
    print(np.min(significance))

    ax = axes[i]
    ax.set_global()
    ax.coastlines(linewidth=0.3)

    contour = ax.contourf(corr_diff['lon'], corr_diff['lat'], corr_diff,
                            levels=bounds, cmap=cmap, norm=norm, transform=data_proj)
                            
    result = significance < 0.05
    print(result)
                            
    ax.contourf(result['lon'], result['lat'], np.ma.masked_where(~result, result), 
        levels=[0,1], colors='none', hatches=['...'], transform=data_proj)
    
    ax.text(0.02, 0.96, fig_labels[i], transform=ax.transAxes,
        fontsize=10, fontweight='bold', va='top', ha='left',
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))
    
           
# --- Map Styling ---
apply_multipanel_map_styling(fig, axes_3d, contour, bounds, label_font_size=10, colorbar_fraction=0.03)


# --- Saving output ---
output_label = f"/unity/f1/cemerson/CODE/interp_{target_enso}_comparison_test"
fig.savefig(output_label + '.png', dpi=600, bbox_inches="tight") 