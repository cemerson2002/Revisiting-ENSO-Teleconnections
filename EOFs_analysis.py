# produces multi-panel correlation figures for all identified temporal oscilations
# (i.e. synoptic, subseasonal, interannual, decadal, and secular)

from supp_functions.preproc_functions import identify_wavelengths, read_enso_data
from supp_functions.calc_functions import corr_significance, format_data, BH_test
from supp_functions.plotting_functions import apply_multipanel_map_styling

from matplotlib.colors import ListedColormap
import matplotlib
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import regionmask
import os
import json

# --- INITIAL SYSTEM CONFIGURATION ---
directory = "/unity/f1/cemerson/data/"
filter_file = directory + "fil_sfc_temp/filtered_CRU_temp.nc"
BH = True
enso_index = ['ONI', 'MEI', 'SOI']

variable = 'temp'
label = '/unity/f1/cemerson/'+variable+'_BH_EOFs'
dataset = filter_file.split('_')[3]

colors = ['navy', 'royalblue', 'lightblue', 'white', 'lightcoral', 'crimson', 'darkred']
bounds = [-1, -0.8, -0.5, -0.2, 0.2, 0.5, 0.8, 1]
cmap = ListedColormap(colors)

fig_proj = ccrs.PlateCarree()
data_proj = ccrs.PlateCarree()

matplotlib.rcParams['axes.linewidth'] = 0.5
land = regionmask.defined_regions.natural_earth_v5_0_0.land_110
xr.set_options(use_flox=True)


# --- PRE-PROCESSING ---
dataset_ds = xr.open_dataset(filter_file)
dataset_da = dataset_ds[list(dataset_ds.data_vars)[0]]

# --- WAVELENGTH ANALYSIS ---
print('Starting wavelength analysis tracking...')
 
modes_file = directory + variable+"_"+dataset+"_modes.json"
print(modes_file)
 
if os.path.exists(modes_file):
    print("Found saved data. Loading from file...")
    with open(modes_file, 'r') as f:
        mode_dict = json.load(f)

else:
    print("No saved data found. Running function...")
 
    mode_dict = identify_wavelengths(dataset_da)
 
    with open(modes_file, 'w') as f:
        json.dump(mode_dict, f, indent=4)

print(mode_dict)

print("Executing one-time global data array seasonal weighting...")
processed_full_dataset = format_data(dataset_da)
print(processed_full_dataset)

enso_data_dict = {index: read_enso_data(index, interannual=True) for index in enso_index}

# --- DATA COLLECTION ---

temporal_definitions = [
    {"label": "Only Interannual", "exclude": None, "only": "Interannual"},
#    {"label": "Synoptic removed", "exclude": ["Synoptic"], "only": None},
#    {"label": "Subseasonal removed", "exclude": ["Subseasonal"], "only": None},
    {"label": "Decadal removed", "exclude": ["Decadal"], "only": None},
    {"label": "Secular removed", "exclude": ["Secular"], "only": None}
]

col, row = len(enso_index), len(temporal_definitions)
fig_labels = [chr(i) for i in range(ord('a'), ord('a') + row * col)]

fig, axes_3d = plt.subplots(ncols=col, nrows=row, subplot_kw={'projection': fig_proj},
                            figsize=(col*2.6, row*1.6), constrained_layout=True, 
                            sharex='col', sharey='row')


for i, period in enumerate(temporal_definitions):
    p_label = period["label"]
    print(f"Calculating: {p_label}")
    
    if period["only"] is not None:
        active_modes = mode_dict.get(period["only"], [])
    else:
        active_modes = []
        for key, modes in mode_dict.items():
            if not any(excl in key for excl in period["exclude"]):
                active_modes.extend(modes if isinstance(modes, list) else [modes])

    print(active_modes)
    sub_dataset = processed_full_dataset.isel(modes=active_modes).sum(dim='modes')
    
    for j, index in enumerate(enso_index):
        print(index)
        print(i)
        print(j)

        ax = axes_3d[i,j]
        
        ax.set_global()
        ax.coastlines(linewidth=0.3)
        
        panel_idx = i * col + j
        print(panel_idx)
        
        ax.text(0.02, 0.96, fig_labels[panel_idx], transform=ax.transAxes,
                    fontsize=9, fontweight='bold', va='top', ha='left',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))
    
        corr1, pvals, _ = corr_significance(sub_dataset, enso_data_dict[index])
        
        pval_BH = BH_test(pvals)
        sig_BH = pval_BH < 0.05
        
        if variable in ['temp', 'precip']:
            mask = land.mask(sub_dataset['lon'], sub_dataset['lat']) 
            
            corr1 = corr1.where(mask == 0)
            sig_BH = sig_BH.where(mask == 0).fillna(False).astype(bool)
 
        contour = ax.contourf(corr1['lon'], corr1['lat'], corr1,
                              levels=bounds, cmap=cmap, transform=data_proj)
 
        
        ax.contourf(sig_BH['lon'], sig_BH['lat'], np.ma.masked_where(~sig_BH, sig_BH), 
                    levels=[0, 1], colors='none', hatches=['...'], transform=data_proj)

apply_multipanel_map_styling(fig=fig, axes_3d=axes_3d, col_labels=enso_index, 
    contour=contour, bounds=bounds, colorbar_fraction = 0.03, label_font_size=6)

fig.savefig(label + '.png', dpi=600, bbox_inches="tight")