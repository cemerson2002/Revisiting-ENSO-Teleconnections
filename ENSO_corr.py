import xarray as xr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t

from supp_functions.preproc_functions import read_enso_data

def process_DJF_seasonal_mean(da_or_ds):
    """
    Groups monthly data into DJF seasons (Dec-Jan-Feb) and calculates the 
    weighted seasonal average based on the number of days in each month.
    December of the prior year is grouped into the current year's winter.
    """
    # 1. Isolate December, January, and February
    ds_djf = da_or_ds.sel(time=da_or_ds['time'].dt.month.isin([12, 1, 2]))
    
    # 2. Shift December forward by 1 month so that Dec 2020, Jan 2021, and Feb 2021 
    # all share the same 'year' coordinate (2021)
    shifted_year = (ds_djf['time.year'] + (ds_djf['time.month'] == 12).astype(int)).rename('year')
    ds_djf.coords['year'] = shifted_year

    # 3. Apply time-weighted grouping across the 'year' coordinate dimension
    weighted_djf = ds_djf.groupby('year').apply(
        lambda group: group.weighted(group.time.dt.days_in_month).mean(dim='time')
    )
    return weighted_djf


filtered_file = "/unity/f1/cemerson/data/ENSO_filter/ENSO_Interannual.csv"

labels = ['a.', 'b.']

enso_index = ['ONI', 'MEI', 'SOI_other_raw']

ENSO_list = []

alpha = 0.05

for index in enso_index:
    ENSO_list.append(read_enso_data(index))
unfiltered = xr.merge(ENSO_list)
unfiltered_df = unfiltered.to_dataframe()
unfiltered_df = unfiltered_df.rename(columns={'SOI_other_raw':'SOI'})

data = pd.read_csv(filtered_file, index_col='time')
data = data.to_xarray()
data['time'] = pd.to_datetime(data['time'].values)
filtered = process_DJF_seasonal_mean(data)
filtered_df = filtered.to_dataframe()
filtered_df = filtered_df[['ONI','MEI','SOI']]

print(unfiltered_df)
print(filtered_df)

# select only overlapping years (1979-2005)
filtered_df = filtered_df.loc['1979':'2024']
unfiltered_df = unfiltered_df.loc['1979':'2024']

print(unfiltered_df)
print(filtered_df)

ENSO = [unfiltered_df, filtered_df]

fig, axes = plt.subplots(1,2, figsize=(12,5)) #, constrained_layout=True)
mappable = None

for label, index, ax in zip(labels, ENSO, axes):
    ax = ax
    index = index.dropna()
    corr = index.corr()
    n = index.shape[0]

    t_stat = corr * np.sqrt((n-2) / (1-corr**2))
    pvals = 2*(1-t.cdf(np.abs(t_stat), df=n-2))

    pval_df = pd.DataFrame(pvals, index=index.columns, columns=index.columns)

    annot = corr.copy().astype(str)

    for i in corr.index:
        for j in corr.columns:
            val = corr.loc[i,j]
            if val != 1:
                if pval_df.loc[i,j] < alpha:
                    annot.loc[i,j] = f"{val:.2f}*"
                else:
                    annot.loc[i,j] = f"{val:.2f}"
            else:
                annot.loc[i,j] = f"{val:.2f}"

    sns.heatmap(corr, ax=ax, annot=annot, fmt="", cmap='coolwarm', 
            vmin=-1, vmax=1, center=0, cbar=False)
    ax.text(-0.5,-0.1,label, weight='bold')
    #sns.pairplot(df)

    ax.tick_params(axis='both', which='both', length=0)

    if mappable is None:
        mappable = ax.collections[0]


fig.colorbar(mappable, ax=axes, orientation='vertical', 
        fraction = 0.05, pad=0.02)

plt.savefig('ENSOcorrs.png', dpi=600)
