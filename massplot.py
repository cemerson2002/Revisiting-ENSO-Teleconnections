# produces multi-panel figures of ENSO teleconnection correlations from "Mode 1" of EEMD

from supp_functions.preproc_functions import read_enso_data, format_data
from supp_functions.calc_functions import corr_significance, BH_test
from supp_functions.plotting_functions import apply_multipanel_map_styling

import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.mpl.ticker as cticker
import matplotlib.ticker as mticker
import regionmask

# --- GLOBAL MATPLOTLIB STYLING ---
matplotlib.rcParams['axes.linewidth'] = 0.5

land = regionmask.defined_regions.natural_earth_v5_0_0.land_110

#reanal_directory = "/unity/f1/cemerson/data/200_filtereed/"
#reanal_directory = "/unity/f1/cemerson/data/500_filtered/"
#reanal_directory = "/unity/f1/cemerson/data/filtered_precip/"
reanal_directory = "/unity/f1/cemerson/data/temp_filtered/"
#reanal_directory = "/unity/f1/cemerson/data/surface_filtered/"
#reanal_directory = "/unity/f1/cemerson/data/fil_sfc_temp/"
reanal_files = sorted(os.listdir(reanal_directory))

enso_index = ['ONI', 'MEI', 'SOI']
LABELS = ["ERA5", "JRA55", "MERRA2", "NCEP-R2"]
#LABELS = ["CHIRPS", "CRU", "3 mo. lag\nCHIRPS", "3 mo. lag\nCRU"]
#LABELS = ["CRU", "3 mo. lag\nCRU"]

orig_var = 'TEMP'
met_is_filter = ['All']
BH = True

colors = ['navy', 'royalblue', 'lightblue', 'white', 'lightcoral', 'crimson', 'darkred']
bounds = [-1, -0.8, -0.5, -0.2, 0.2, 0.5, 0.8, 1]

cmap = ListedColormap(colors)
norm = BoundaryNorm(bounds, cmap.N)

# Set dynamic output string name based on metadata
var_suffix = f"{orig_var}_BH" if BH else f"{orig_var}_ttest"
if met_is_filter:
    var_suffix += "_metIMF"
label = os.path.join('/unity/f1/cemerson/CODE/', var_suffix)

fig_proj = ccrs.PlateCarree()
data_proj = ccrs.PlateCarree()

row, col = len(reanal_files), len(enso_index)
fig_labels = [chr(i) for i in range(ord('a'), ord('a') + row * col)]

fig, axes_3d = plt.subplots(ncols=col, nrows=row, subplot_kw={'projection': fig_proj},
                            figsize=(col * 2.6, row * 1.2), constrained_layout=True, 
                            sharex='col', sharey='row')

lags = ['DJF'] 

# (file_idx, lag_idx) -> row_index
#row_map = {
#    (0, 0): 0,  # CHIRPS regular
#    (1, 0): 1,  # CRU regular
#    (0, 1): 2,  # CHIRPS lag
#    (1, 1): 3   # CRU lag
#}

for i, filename in enumerate(reanal_files):
    print(f'Processing file {i}: {filename}')
    reanal_filepath = os.path.join(reanal_directory, filename)
    
    with xr.open_dataset(reanal_filepath) as ds:
        var_name = list(ds.data_vars)[0]
        dataset_ds = ds[var_name][0, :, :, :]
        
    processed_ds = format_data(dataset_ds)
 
    mask = None
    if orig_var in ['TEMP', 'PRECIP']:
        mask = land.mask(processed_ds['lon'], processed_ds['lat'])
        
    for lag_idx, lag in enumerate(lags):
            
        # MATH MAGIC: This strictly forces:
        # lag_idx=0, i=0 -> Row 0 (CHIRPS regular)
        # lag_idx=0, i=1 -> Row 1 (CRU regular)
        # lag_idx=1, i=0 -> Row 2 (CHIRPS 3-month lag)
        # lag_idx=1, i=1 -> Row 3 (CRU 3-month lag)
        #row = row_map[(i, lag_idx)]
        
        for j, index in enumerate(enso_index):
            
            
            ENSO_da = read_enso_data(index, interannual=True, months=lag) 
            
            ax = axes_3d[i, j]
            ax.set_global()
            ax.coastlines(linewidth=0.3)
            
            panel_idx = i * len(enso_index) + j
            ax.text(0.02, 0.96, fig_labels[panel_idx], transform=ax.transAxes,
                    fontsize=8, fontweight='bold', va='top', ha='left',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))
            
            # Run correlations and plot
            corr1, pvals, _ = corr_significance(processed_ds, ENSO_da)
            pval_BH = BH_test(pvals)
            sig_BH = pval_BH < 0.05
            
            if orig_var in ['TEMP', 'PRECIP']:
                corr1 = corr1.where(mask == 0)
                sig_BH = sig_BH.where(mask == 0).fillna(False).astype(bool)
     
            contour = ax.contourf(corr1['lon'], corr1['lat'], corr1,
                                  levels=bounds, cmap=cmap, transform=data_proj)
     
            
            ax.contourf(sig_BH['lon'], sig_BH['lat'], np.ma.masked_where(~sig_BH, sig_BH), 
                        levels=[0, 1], colors='none', hatches=['...'], transform=data_proj)
            
# Run the reusable formatting engine
apply_multipanel_map_styling(
    fig=fig, 
    axes_3d=axes_3d, 
    row_labels=LABELS, 
    col_labels=enso_index, 
    contour=contour, 
    bounds=bounds,
    label_font_size=6, 
    colorbar_fraction=0.03
)

fig.savefig(label + '.png', dpi=600, bbox_inches="tight")
