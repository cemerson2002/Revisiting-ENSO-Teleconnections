# produces bar chart figures displaying change in significant correlation area

from supp_functions.preproc_functions import identify_wavelengths, read_enso_data
from supp_functions.calc_functions import corr_significance, format_data, BH_test

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import regionmask
import pandas as pd
import matplotlib.patches as mpatches
import sys
import os

import json


# --- INITIAL SYSTEM CONFIGURATION ---
#filter_file = "/unity/f1/cemerson/data/200_filtereed/ERA5_filtered_height_200.nc"
filter_file = "/unity/f1/cemerson/data/surface_filtered/filtered_CHIRPS_precip.nc"
#filter_file = "/unity/f1/cemerson/data/fil_sfc_temp/filtered_CRU_temp.nc"
BH = True
enso_index = ['ONI', 'MEI', 'SOI']
variable = 'PRECIP'
dataset = 'CHIRPS'

modes_file = '/unity/f1/cemerson/data/'+variable+'_'+dataset+'_modes.json'


def compute_bin_percentages(da, enso_data_dict, enso_index, region_filters, land_mask, use_BH=True):
    corr_storage = {}
    for index in enso_index:
        print(index)
        corr, pvals, _ = corr_significance(da, enso_data_dict[index])

        pval_BH = BH_test(pvals)
        sig = pval_BH < 0.05
        
        if index == 'SOI':
            corr = -1 * corr
            
        corr_masked = corr.where(sig, 0.0).where(~np.isnan(corr))
        corr_storage[index] = corr_masked
        

    results = {}
    for reg_name, filter_func in region_filters.items():
        for index in enso_index:
            reg_vals = filter_func(corr_storage[index], land_mask).values
            reg_vals = reg_vals[~np.isnan(reg_vals)]
            total    = len(reg_vals) if len(reg_vals) > 0 else 1
            results[(reg_name, index)] = {
                'St_Neg': np.sum(reg_vals <  -0.5)                          / total * 100,
                'Wk_Neg': np.sum((reg_vals >= -0.5) & (reg_vals <=  -0.2))  / total * 100,
                'Neu':    np.sum((reg_vals > -0.2) & (reg_vals <   0.2))   / total * 100,
                'Wk_Pos': np.sum((reg_vals >=  0.2) & (reg_vals <=   0.5))  / total * 100,
                'St_Pos': np.sum(reg_vals >  0.5)                          / total * 100,
            }

    return results
    
def get_active_modes(mode_dict, exclude=None, only=None):
    always_exclude = {'Original', 'Noise'}
    if only is not None:
        return list(mode_dict.get(only, []))
    active = []
    for key, modes in mode_dict.items():
        if key in always_exclude:
            continue
        if exclude and any(excl in key for excl in exclude):
            continue
        active.extend(modes if isinstance(modes, list) else [modes])
    return active

# --- PRE-PROCESSING ---
dataset_ds = xr.open_dataset(filter_file)
dataset_da = dataset_ds[list(dataset_ds.data_vars)[0]]

print('Starting wavelength analysis tracking...')
# Check if the file already exists
if os.path.exists(modes_file):
    print("Found saved data. Loading from file...")
    with open(modes_file, 'r') as f:
        mode_dict = json.load(f)

else:
    print("No saved data found. Running function...")
    
    mode_dict = identify_wavelengths(dataset_da)
    
    # Save the output to a JSON file for next time
    with open(modes_file, 'w') as f:
        json.dump(mode_dict, f, indent=4)
    print("Data saved successfully for next time!")

print(mode_dict)


print("Executing one-time global data array seasonal weighting...")
processed_full_dataset = format_data(dataset_da)
land_mask = regionmask.defined_regions.natural_earth_v5_0_0.land_110.mask(processed_full_dataset['lon'], processed_full_dataset['lat'])
enso_data_dict = {index: read_enso_data(index, interannual=True) for index in enso_index}

# ==========================================
# REGION DEFINITIONS (using existing land_mask)
# ==========================================


region_filters = {
    'Global':            lambda c, m: c.where(m.values == 0), # only land for precip and temp
    #'Trop. Ocean':    lambda c, m: c.where(m.values != 0).where(np.abs(c.lat) <= 23.5),
    'Trop. Land':     lambda c, m: c.where(m.values == 0).where(np.abs(c.lat) <= 23.5),
    #'NH Midlat Ocean':  lambda c, m: c.where(m.values != 0).where((c.lat > 23.5)  & (c.lat <= 60)),
    'NH Midlat Land':   lambda c, m: c.where(m.values == 0).where((c.lat > 23.5)  & (c.lat <= 60)),
    #'NH Polar Ocean': lambda c, m: c.where(m.values != 0).where(c.lat > 60),
    #'NH Polar Land':  lambda c, m: c.where(m.values == 0).where(c.lat > 60),
    #'SH Midlat Ocean':  lambda c, m: c.where(m.values != 0).where((c.lat < -23.5) & (c.lat >= -60)),
    'SH Midlat Land':   lambda c, m: c.where(m.values == 0).where((c.lat < -23.5) & (c.lat >= -60)),
    #'SH Polar Ocean': lambda c, m: c.where(m.values != 0).where(c.lat < -60),
    #'SH Polar Land':  lambda c, m: c.where(m.values == 0).where(c.lat < -60),
}

# Preserve lat coordinate through filter
for k in region_filters:
    region_filters[k] = lambda c, m, f=region_filters[k]: f(c, m).assign_coords(lat=c.lat)

# ==========================================
# DATA COLLECTION
# ==========================================

# --- MODE 0 BASELINE ---
baseline_da = processed_full_dataset.isel(modes=0)

# --- TEMPORAL PERIOD DEFINITIONS ---

temporal_definitions = [
    {"label": "Interannual Only",      "only": "Interannual", "exclude": None},
    {"label": "Synoptic Removed",      "only": None, "exclude": ["Synoptic"]},
    {"label": "Subseasonal Removed",   "only": None, "exclude": ["Subseasonal"]},
    {"label": "Decadal Removed",       "only": None, "exclude": ["Decadal"]},
    #{"label": "Secular Removed",       "only": None, "exclude": ["Secular"]},
]

bin_cols = ['St_Neg', 'Wk_Neg', 'Neu', 'Wk_Pos', 'St_Pos']

print("Computing mode 0 baseline correlations...")
baseline_pcts = compute_bin_percentages(
    baseline_da, enso_data_dict, enso_index, region_filters, land_mask, use_BH=BH
)
print(baseline_pcts)

print("Computing temporal period deltas...")
delta_records = []

for period in temporal_definitions:
    p_label = period["label"]
    print(f"  Computing: {p_label}")
    active_modes = get_active_modes(mode_dict, exclude=period["exclude"], only=period["only"])

    if not active_modes:
        print(f"  Skipping {p_label}: no active modes found")
        continue

    try:
        period_da   = processed_full_dataset.isel(modes=active_modes).sum(dim='modes')
        period_pcts = compute_bin_percentages(
            period_da, enso_data_dict, enso_index, region_filters, land_mask, use_BH=BH
        )
        for (reg_name, index), bins in period_pcts.items():
            base = baseline_pcts[(reg_name, index)]
            record = {'Period': p_label, 'Region': reg_name, 'Index': index}
            for b in bin_cols:
                record[b] = bins[b] - base[b]
            delta_records.append(record)

    except Exception as e:
        print(f"  Failed {p_label}: {e}")

df_delta = pd.DataFrame(delta_records)
print(df_delta)

# ==========================================
# STACKED BAR CHART (OMITTING NEUTRAL)
# ==========================================

delta_period_order = [p["label"] for p in temporal_definitions]

#row_definitions = [
#    ('Global',     'Global',          None),
#    ('NH Polar',   'NH Polar Land',   'NH Polar Ocean'),
#    ('NH Midlat',  'NH Midlat Land',  'NH Midlat Ocean'),
#    ('Tropics',   'Trop. Land',      'Trop. Ocean'),
#    ('SH Midlat',  'SH Midlat Land',  'SH Midlat Ocean'),
#    ('SH Polar',   'SH Polar Land',   'SH Polar Ocean'),
#]

row_definitions = [
    ('Global',     'Global',          None),
    #('NH Polar',   'NH Polar Land',   None),
    ('NH Midlat',  'NH Midlat Land',  None),
    ('Tropics',   'Trop. Land',      None),
    ('SH Midlat',  'SH Midlat Land',  None)
    #('SH Polar',   'SH Polar Land',   None)
]

bar_colors = {
    'St_Neg': '#1a3a6b',
    'Wk_Neg': '#5a90d9',
    'Wk_Pos': '#d96a5a',
    'St_Pos': '#6b1a1a',
}

categories  = ['St_Neg', 'Wk_Neg', 'Wk_Pos', 'St_Pos']
n_rows      = len(row_definitions)
n_cols      = len(enso_index)
n_periods   = len(delta_period_order)
 
# Bar geometry — two bars per period group, with a small gap between them
BAR_W   = 0.38   # width of each individual bar
BAR_GAP = 0.05   # gap between the Land and Ocean bar within a group
# Group centres sit at integer x positions
x_centres    = np.arange(n_periods, dtype=float)
land_offsets  = x_centres - (BAR_W / 2 + BAR_GAP / 2)
ocean_offsets = x_centres + (BAR_W / 2 + BAR_GAP / 2)
 
fig, axes = plt.subplots(
    nrows=n_rows, ncols=n_cols,
    figsize=(4.2 * n_cols, 2 * n_rows),
    sharex=True, sharey=True,
    constrained_layout=True,
)
 
fig_labels = [chr(i) for i in range(ord('a'), ord('a') + n_rows * n_cols)]
 
 
def draw_stacked_bars(ax, region, x_offsets, df_delta, hatch=None, alpha=1.0):
    """Draw one set of stacked diverging bars for *region* at *x_offsets*."""
    if region is None:
        return
    subset = (
        df_delta[(df_delta['Region'] == region) & (df_delta['Index'] == idx)]
        .set_index('Period')
        .reindex(delta_period_order)
        .fillna(0)
    )
    pos_base = np.zeros(n_periods)
    neg_base = np.zeros(n_periods)
    kw = dict(width=BAR_W, linewidth=0.4, alpha=alpha,
              hatch=hatch, edgecolor='white' if hatch is None else '#555555')
 
    for cat in categories:
        vals = subset[cat].values
        bottoms = np.zeros(n_periods)
        for p in range(n_periods):
            val = vals[p]
            if val >= 0:
                bottoms[p] = pos_base[p]
                pos_base[p] += val
            else:
                bottoms[p] = neg_base[p]
                neg_base[p] += val
        ax.bar(x_offsets, vals, bottom=bottoms, color=bar_colors[cat], **kw)
 
 
for r, (row_label, land_region, ocean_region) in enumerate(row_definitions):
    paired = ocean_region is not None   # True when we have both Land and Ocean
 
    for c, idx in enumerate(enso_index):
        ax = axes[r][c]
 
        if paired:
            draw_stacked_bars(ax, land_region,  land_offsets,  df_delta, hatch=None,  alpha=1.0)
            draw_stacked_bars(ax, ocean_region, ocean_offsets, df_delta, hatch='//', alpha=0.85)
        else:
            # Single region (Global) — centre bars on x_centres
            draw_stacked_bars(ax, land_region, x_centres, df_delta, hatch=None, alpha=1.0)
 
        # Zero line and grid
        ax.axhline(0, color='black', linewidth=1.0, zorder=5)
        ax.xaxis.grid(True, linestyle=':', alpha=0.3, color='gray', zorder=0)
        ax.yaxis.grid(True, linestyle='--', alpha=0.3, color='gray', zorder=0)
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
 
        # Column header (ENSO index)
        if r == 0:
            ax.set_title(idx, fontweight='bold', fontsize=15, pad=10)
 
        # Row label on the right
        if c == n_cols - 1:
            ax.text(1.04, 0.5, row_label, transform=ax.transAxes,
                    rotation=270, ha='left', va='center',
                    fontsize=13, fontweight='bold', color='#333333')
 
        # x-axis tick labels only on the bottom row
        ax.set_xticks(x_centres)
        if r == n_rows - 1:
            ax.set_xticklabels(delta_period_order, fontsize=10, rotation=25, ha='right')
            ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
        else:
            ax.tick_params(axis='x', which='both', bottom=True, labelbottom=False)
 
        if c == 0:
            ax.set_ylabel('Change in Area (%)', fontsize=10)
 
        ax.set_ylim(-55, 55)
        ax.set_yticks([-50, -25, 0, 25, 50])
        ax.tick_params(axis='y', labelsize=10)
 
        # Panel label (a, b, c …)
        panel_idx = r * n_cols + c
        ax.text(0.02, 0.96, fig_labels[panel_idx], transform=ax.transAxes,
                fontsize=13, fontweight='bold', va='top', ha='left',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=2))
 
# ---- Legend ----
# Colour legend (correlation bins)
color_patches = [
    mpatches.Patch(color=bar_colors['St_Neg'], label='Strong Negative  (r = -0.5)'),
    mpatches.Patch(color=bar_colors['Wk_Neg'], label='Weak Negative  (-0.5 < r < -0.2)'),
    mpatches.Patch(color=bar_colors['Wk_Pos'], label='Weak Positive  (0.2 = r < 0.5)'),
    mpatches.Patch(color=bar_colors['St_Pos'], label='Strong Positive  (r = 0.5)'),
]
 
# Hatch legend (surface type)
land_patch  = mpatches.Patch(facecolor='#aaaaaa', edgecolor='white',   linewidth=0.4,
                              label='Land (solid)')
ocean_patch = mpatches.Patch(facecolor='#aaaaaa', edgecolor='#555555', linewidth=0.4,
                              hatch='//', alpha=0.85, label='Ocean (hatched)')
 
all_handles = color_patches + [
    mpatches.Patch(visible=False),   # spacer
    land_patch, ocean_patch,
]
 
fig.legend(
    handles=all_handles,
    loc='lower center', ncol=6,
    bbox_to_anchor=(0.5, -0.05),
    fontsize=10, frameon=False,
)
 
fig.savefig(variable+'_delta_diverging_bars_by_region.png', dpi=600, bbox_inches='tight')
plt.show()