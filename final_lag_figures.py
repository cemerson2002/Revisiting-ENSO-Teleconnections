"""
Multi-Reanalysis ENSO Correlation Plotter
=========================================
Generates an 8-panel plot (4 rows x 2 cols) comparing 4 reanalysis datasets 
(ERA5, JRA55, MERRA2, NCEPR2) against two different ENSO indices.
The left column uses the DJF ENSO index, and the right column uses an easily 
configurable ENSO index. The target climate variable is fixed to the DJF season.
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
from matplotlib.colors import ListedColormap, BoundaryNorm
from cartopy.util import add_cyclic_point

# Import from your custom modules
from supp_functions.calc_functions import format_data, corr_significance
from supp_functions.plotting_functions import apply_multipanel_map_styling

warnings.filterwarnings("ignore", category=RuntimeWarning)

# -----------------------------------------------------------------------------
# 1. CONFIGURATION & FILE PATHS
# -----------------------------------------------------------------------------

# --- EASY TOGGLE FOR COLUMNS AND LAGS ---
LEFT_COL_SEASON  = 'DJF' 
LEFT_LAG_YEARS   = 0     # 0 = Same climate year

RIGHT_COL_SEASON = 'SON' # <-- Set to SON
RIGHT_LAG_YEARS  = 1     # 1 = Shift forward by 1 year (SON of previous year aligns with DJF)

# Dictionary mapping reanalysis names to their respective file paths.
DATASETS = {
    'CRU_1':   "/unity/f1/cemerson/data/surface_filtered/filtered_CRU_precip.nc",
    'CHIRPS':  "/unity/f1/cemerson/data/surface_filtered/filtered_CHIRPS_precip.nc", 
    'CRU_2': "/unity/f1/cemerson/data/fil_sfc_temp/filtered_CRU_temp.nc"
}

ENSO_EXCEL_PATH = "/unity/f1/cemerson/data/ENSO_indexes.xlsx"
TARGET_SHEETS    = ["ONI", "MEI", "SOI"]
ALPHA           = 0.05  # Significance level for BH-FDR

# -----------------------------------------------------------------------------
# 2. DATA PREPARATION
# -----------------------------------------------------------------------------

print("Loading and formatting Reanalysis datasets (DJF)...")
clm_djf_dict = {}

for name, path in DATASETS.items():
    print(f" -> Processing {name}...")
    ds = xr.open_dataset(path)
    
    # Dynamically grab the main data variable
    var_name = list(ds.data_vars)[0]
    da = ds[var_name][0,:,:,:]
    
    # 1. Rigorously standardize coordinate names 
    # (Accounts for 'Time', 'latitude', 'longitude', etc.)
    rename_dict = {}
    for coord in da.coords:
        if 'time' in str(coord).lower() and coord != 'time':
            rename_dict[coord] = 'time'
        elif 'lat' in str(coord).lower() and coord != 'lat':
            rename_dict[coord] = 'lat'
        elif 'lon' in str(coord).lower() and coord != 'lon':
            rename_dict[coord] = 'lon'
    if rename_dict:
        da = da.rename(rename_dict)
        
    # Format directly to DJF using your robust format_data function
    clm_djf_dict[name] = format_data(da, season='DJF')

# -----------------------------------------------------------------------------
# 3. ENSO INDEX LOOP & PLOTTING
# -----------------------------------------------------------------------------
# We now loop over each ENSO index, load its data, and immediately create its plot.
for enso_idx in TARGET_SHEETS:
    print(f"\nProcessing ENSO Index: {enso_idx} ({LEFT_COL_SEASON} and {RIGHT_COL_SEASON})...")
    
    # Read ENSO Data for the current index
    df_enso = pd.read_excel(ENSO_EXCEL_PATH, sheet_name=enso_idx).set_index('year')

    # Extract the fixed left column and dynamic right column
    enso_left  = df_enso[LEFT_COL_SEASON].dropna()
    enso_right = df_enso[RIGHT_COL_SEASON].dropna()
    
    # Ensure indices are integer type
    enso_left.index  = enso_left.index.astype(int)
    enso_right.index = enso_right.index.astype(int)
    
      # --- DYNAMIC LAG ALIGNMENT ---
    enso_left.index  = enso_left.index + LEFT_LAG_YEARS
    enso_right.index = enso_right.index + RIGHT_LAG_YEARS
    
    # Convert pandas Series to xarray DataArrays
    enso_left_da  = xr.DataArray(enso_left.values, coords=[enso_left.index], dims=['year'])
    enso_right_da = xr.DataArray(enso_right.values, coords=[enso_right.index], dims=['year'])

    print(f"Calculating correlations and building 8-panel plot for {enso_idx}...")

    # Set up the 4x2 grid
    fig, axes = plt.subplots(
        nrows=3, ncols=2, 
        figsize=(9, 6), 
        subplot_kw={'projection': ccrs.PlateCarree()}, 
        constrained_layout=True 
    )

    # Set labels dynamically based on the current index and configured seasons
    col_labels = [
        f'{enso_idx} ({LEFT_COL_SEASON}) vs Climate (DJF)', 
        f'{enso_idx} ({RIGHT_COL_SEASON}) vs Climate (DJF)'
    ]
    row_labels = list(DATASETS.keys())

    # The target order maps directly to the columns [Left (Col 0), Right (Col 1)]
    enso_targets = [enso_left_da, enso_right_da]

    # Plotting parameters
    contour_obj = None
    fig_labels = [chr(i) for i in range(ord('a'), ord('a') + 8)]
    colors = ['navy', 'royalblue', 'lightblue', 'white', 'lightcoral', 'crimson', 'darkred']
    bounds = [-1, -0.8, -0.5, -0.2, 0.2, 0.5, 0.8, 1]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)

    for i, (name, clm_da) in enumerate(clm_djf_dict.items()):
        for j, enso_da in enumerate(enso_targets):
            ax = axes[i, j]
            
            # 1. Calculate correlation and BH-corrected significance
            corr, sig_mask, pval, sig_bh = corr_significance(clm_da, enso_da, alpha=ALPHA)

            lons = corr.lon.values
            lats = corr.lat.values
            
            # Add cyclic point to the correlation data
            corr_cyclic, lons_cyclic = add_cyclic_point(corr.values, coord=lons)
            
            # 2. Draw map contour
            contour_obj = ax.contourf(
                lons_cyclic, lats, corr_cyclic, 
                transform=ccrs.PlateCarree(),
                levels=bounds, cmap=cmap, extend='both'
            )
            
            # 3. Add Cartopy features
            ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor='#333333')
            
            panel_idx = i * 2 + j
            ax.text(0.02, 0.96, fig_labels[panel_idx], transform=ax.transAxes,
                    fontsize=9, fontweight='bold', va='top', ha='left',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))
            
            # Add cyclic point to the significance mask (prevents white line gap)
            sig_cyclic, sig_lons_cyclic = add_cyclic_point(sig_bh.astype(int).values, coord=lons)
            
            # 4. Significance Stippling (using the cyclic BH corrected mask)
            ax.contourf(sig_lons_cyclic, lats, sig_cyclic, 
                        levels=[0.5, 1.5], colors='none', hatches=['...'], transform=ccrs.PlateCarree())
                
            # Ensure extent is global or matched to data
            ax.set_global()

    # 5. Apply custom map styling
    apply_multipanel_map_styling(
        fig=fig, 
        axes_3d=axes, 
        col_labels=col_labels, 
        contour=contour_obj, 
        bounds=bounds, 
        row_labels=row_labels, 
        label_font_size=9,
        colorbar_fraction=0.03
    )

    # Save specific to the current ENSO index
    output_filename = f"enso_lagprecip_{enso_idx}_{LEFT_COL_SEASON}_vs_{RIGHT_COL_SEASON}.png"
    plt.savefig(output_filename, dpi=200)
    print(f"Plot saved successfully to {output_filename}")
    
    # Close figure to free memory before next iteration
    plt.close(fig)