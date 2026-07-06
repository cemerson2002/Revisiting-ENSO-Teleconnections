# produces PC sensitivity figures

from supp_functions.preproc_functions import read_enso_data, format_data, identify_wavelengths
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
import regionmask
import json

from dask.diagnostics import ProgressBar

# --- GLOBAL MATPLOTLIB STYLING ---
matplotlib.rcParams['axes.linewidth'] = 0.5
land = regionmask.defined_regions.natural_earth_v5_0_0.land_110
xr.set_options(use_flox=True)

fig_proj = ccrs.PlateCarree()
data_proj = ccrs.PlateCarree()

colors = ['navy', 'royalblue', 'lightblue', 'white', 'lightcoral', 'crimson', 'darkred']
bounds = [-1, -0.8, -0.5, -0.2, 0.2, 0.5, 0.8, 1]

cmap = ListedColormap(colors)
norm = BoundaryNorm(bounds, cmap.N)


# --- CONFIGURATIONS ----

directory = "/unity/f1/cemerson/PC_sensitivity/"
datasets = ['CHIRPS', 'CRU']
nPCs = [5, 8, 10, 12]

row, col = len(datasets), len(nPCs)
fig_labels = [chr(i) for i in range(ord('a'), ord('a') + row * col)]

fig, axes_3d = plt.subplots(ncols=col, nrows=row, subplot_kw={'projection': fig_proj},
                            figsize=(col*2.6, row*1.4), constrained_layout=True, 
                            sharex='col', sharey='row')
                   
                            
# --- READ ENSO DATA ---

ENSO_da = read_enso_data('ONI', interannual=True)
                            
  
# --- START PROCESSING LOOP ---
  
for i, dataset in enumerate(datasets):

    for j, count in enumerate(nPCs):
    
        # --- PARSE FILE NAME ---
        try:
            name = dataset.split('_')[0]
            level = dataset.split('_')[1]
        except:
            name = dataset
        
        if name == 'ERA5':
            variable = 'height'
        elif name == 'CHIRPS':
            variable = 'precip'
        elif name == 'CRU':
            variable = 'temp'
            
        file_name = directory + name + '_' + str(count) + '_filtered_' + variable + '.nc'
    
        print(f'Processing file: {file_name}')
        
        
        # --- OPEN DATASET ---
        
        with xr.open_dataset(file_name) as ds:
            var_name = list(ds.data_vars)[0]
            dataset_ds = ds[var_name]
            
        print(dataset_ds['modes'].values)
        
            
        # --- WAVELENGTH ANALYSIS ---
        print('Starting wavelength analysis tracking...')
        
        modes_file = '/unity/f1/cemerson/data/'+variable+'_'+dataset+'_'+str(count)+'_modes.json'
        print(modes_file)
        
        if os.path.exists(modes_file):
            print("Found saved data. Loading from file...")
            with open(modes_file, 'r') as f:
                mode_dict = json.load(f)

        else:
            print("No saved data found. Running function...")
    
            mode_dict = identify_wavelengths(dataset_ds)
    
            with open(modes_file, 'w') as f:
                json.dump(mode_dict, f, indent=4)

        print(mode_dict)
        
        int_modes = mode_dict['Interannual']
        dataset_ds = dataset_ds.isel(modes=int_modes).sum(dim='modes')
        
        with ProgressBar():
            processed_ds = format_data(dataset_ds).load()
        
        print('done formatting')
        # --- BEGIN PLOTTING ---
     
        mask = None
        if variable in ['temp', 'precip']:
            mask = land.mask(processed_ds['lon'], processed_ds['lat']) 
        
        print(i)
        print(j)
        ax = axes_3d[i, j]
        ax.set_global()
        ax.coastlines(linewidth=0.3)
        
        panel_idx = i * len(nPCs) + j 
        print(panel_idx)
        ax.text(0.02, 0.96, fig_labels[panel_idx], transform=ax.transAxes,
                fontsize=8, fontweight='bold', va='top', ha='left',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))
        
        # Run correlations and plot
        print('calculating correlations...')
        corr1, pvals, _ = corr_significance(processed_ds, ENSO_da)
        print(corr1)
        print('conducting BH test...')
        pval_BH = BH_test(pvals)
        sig_BH = pval_BH < 0.05
        
        if variable in ['temp', 'precip']:
            corr1 = corr1.where(mask == 0)
            sig_BH = sig_BH.where(mask == 0).fillna(False).astype(bool)
 
        contour = ax.contourf(corr1['lon'], corr1['lat'], corr1,
                              levels=bounds, cmap=cmap, transform=data_proj)
 
        
        ax.contourf(sig_BH['lon'], sig_BH['lat'], np.ma.masked_where(~sig_BH, sig_BH), 
                    levels=[0, 1], colors='none', hatches=['...'], transform=data_proj)
            
apply_multipanel_map_styling(
    fig=fig, 
    axes_3d=axes_3d, 
    row_labels=datasets, 
    col_labels=nPCs, 
    contour=contour, 
    bounds=bounds,
    label_font_size=6, 
    colorbar_fraction=0.03
)

fig.savefig('/unity/f1/cemerson/PC_other.png', dpi=600, bbox_inches="tight")