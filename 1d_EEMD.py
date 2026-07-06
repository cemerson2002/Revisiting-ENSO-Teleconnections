# conducts one dimensional EEMD on the ENSO index time series to obtain the interannual variations

# --- PACKAGES ---
from PyEMD import EEMD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd

from scipy.signal import find_peaks

# --- CONFIGURATION ---
filename='/unity/f1/cemerson/data/ENSO_indexes.xlsx'
sheet_names = ['ONI', 'JMA_raw', 'MEI_raw', 'MEI_ext_raw', 'SOI_raw']

full_subseasonal = []
full_interannual = []
full_decadal = []
full_secular = []

index_name = ['ONI', 'JMA', 'MEI', 'MEI_ext', 'SOI']

master_dates = pd.date_range(start='1871-01-01', end='2024-12-01', freq='MS')

period_dfs = {
        'Subseasonal': pd.DataFrame(index=master_dates),
        'Interannual': pd.DataFrame(index=master_dates),
        'Decadal': pd.DataFrame(index=master_dates),
        'Secular': pd.DataFrame(index=master_dates)
        }

for name, reg_name in zip(sheet_names, index_name):
    print(name)
    # read in monthly enso index
    df = pd.read_excel(filename, sheet_name=name)

    # convert from 2d to 1d
    df_melted = pd.melt(df, id_vars='year', var_name='month', value_name='index')
    df_melted['date'] = pd.to_datetime(df_melted['year'].astype(str) + 
            '-' + df_melted['month'].astype(str))
    df_melted = df_melted.sort_values('date').reset_index(drop=True)
    df_melted = df_melted[['date', 'index']]
    df_melted = df_melted.dropna()

    # makes array length of monthly indexes for EEMD 
    index = np.array(df_melted['index'], dtype=float)
    years = np.arange(0, len(index), 1)

    # do EEMD
    eemd = EEMD()

    #emd = eemd.EMD
    #emd.extrema_detection='parabol'

    eIMFs = eemd.eemd(index, years)
    eIMFs = eIMFs[2:,:]
    nIMFs = eIMFs.shape[0]
    
    time_periods = {
            'Subseasonal': [],
            'Interannual': [],
            'Decadal': [],
            'Secular': []
            }

    time_periods['Subseasonal'].append(eIMFs[1,:])
    
    for i in range(eIMFs.shape[0]):
        print(f'Mode: {i+2}')
        signal = eIMFs[i]
        signal = pd.Series(signal)  
        peaks, _ = find_peaks(signal, prominence=np.max(signal)*0.005)

        if np.all(signal>=0) or np.all(signal<=0):
            print('secular')
            time_periods['Secular'].append(signal)

        elif len(peaks)<=1:
            print('decadal')
            time_periods['Decadal'].append(signal)

        else:
            min_wavelength = np.percentile(np.diff(peaks), 25)
            max_wavelength = np.percentile(np.diff(peaks), 75)

            print(f'Max: {max_wavelength}')
            print(f'Min: {min_wavelength}')

            if min_wavelength>=1 and max_wavelength<=3: # subseasonal
                print('subseasonal')
                time_periods['Subseasonal'].append(signal)

            elif min_wavelength>3 and max_wavelength<=120: # interannual
                print('interannual')
                time_periods['Interannual'].append(signal)
                
            elif min_wavelength>120: # decadal
                print('decadal')
                time_periods['Decadal'].append(signal)

            else:
                print(f'Mode {i+2} does not fit into any group')
        
    plt.figure(figsize=(12,9))
    ax = plt.subplot(nIMFs+1, 1, 1)
    ax.plot(df_melted['date'], index, 'r')

    for n in range(nIMFs):
        ax = plt.subplot(nIMFs+1, 1, n+2)
        ax.plot(df_melted['date'][:60], eIMFs[n][:60], 'g')
        ax.set_ylabel(f'eIMF {n+1}')

    plt.xlabel('Time')
    plt.tight_layout()
    plt.savefig(f'enso_{reg_name}.png')


    for name, signals in time_periods.items():
        summed = np.sum(signals, axis=0)
        summed_series = pd.Series(summed)

        time_index = df_melted['date'][:len(summed_series)]
        summed_series.index = time_index

        aligned_signal = summed_series.reindex(master_dates)

        period_dfs[name][reg_name] = aligned_signal
        period_dfs[name].index.name = 'time'

for name, df in period_dfs.items():
    df.to_csv(f'ENSO_filter/ENSO_{name}.csv', index=True)

