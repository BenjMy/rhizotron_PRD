#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 15:46:32 2022

@author: ben
"""
import pandas as pd
from matplotlib.pylab import plt
import matplotlib.dates as mdates
import matplotlib.transforms as mtransforms

#%%
imaging_path = '../imaging_ERT_MALM/'
fmt = '%d/%m/%Y,%H:%M'
startD = None
endD = None


#%% Load log files 
# ----------------------------------------------------------------------------


mosaic = '''
            ab
         '''
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=False
                              )


# fig, ax = plt.subplots(1, figsize=(5,3), tight_layout=True)
dleaf = pd.read_csv('../PRD_leaf_parameters - surface sum.csv',
                    # infer_datetime_format=True,
                    thousands=',')


dleaf['datetime'] = pd.to_datetime(dleaf['data'],
                                   format='%d/%m/%Y'
                                   )
dleaf.set_index('datetime',inplace=True)
nameg = []
for g in range(3):
    nameg.append('estimated total LA '+ str(g+1) + ' [cm2]')	


# estimated total LA 1 [cm2]
# estimated total LA 1 [cm2]
legend = ['Model1','Model2','Model3']

markers = ['x', 'o', '^']

dleaf[nameg].plot(
                    # marker= 'v',
                    ylabel='estimated total LA \n (cm2)',
                    label=legend,
                    ax=axs['a'],
                    # cmap='viridis'
                    color = 'k',
                    style=['+-', 'o-', '.--', 's:']   
                    )
axs['a'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
axs['a'].legend(legend)

axs['a'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))


leaf_gsw = pd.read_csv('../leaf_stomatal_conductance_PRD.csv',
                    # infer_datetime_format=True,
                    thousands=',')

leaf_gsw['datetime'] = pd.to_datetime(leaf_gsw['Date'])
leaf_gsw['stress'] = 'low'

leaf_gsw['gsw'] = leaf_gsw['gsw']*1e-3


leaf_gsw['stress'].loc[leaf_gsw['datetime']=='2022-08-06'] = 'high'
# leaf_gsw.groupby(by='datetime').boxplot(column='gsw', ax=axs['b'])
plt.setp(axs['a'].get_xticklabels(), rotation=30, ha='right')


from scipy.stats import ttest_ind

stat, p_value = ttest_ind(leaf_gsw[leaf_gsw['stress']=='high']['gsw'],
                          leaf_gsw[leaf_gsw['stress']=='low']['gsw']
                          )
print(f"t-test: statistic={stat:.4f}, p-value={p_value:.4f}")
# t-test: statistic=-1.5549, p-value=0.1203




import seaborn as sns

# tips = sns.load_dataset("tips")

# Draw a nested boxplot to show bills by day and time
sns.boxplot(x="datetime", y="gsw",
            hue="stress", palette=["m", "g"],
            data=leaf_gsw,
            ax=axs['b'],
            )
axs['b'].set_ylabel(ylabel='gsw \n (mmol m−2 s−1)')
# axs['b'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))

x_dates = leaf_gsw['datetime'].dt.strftime('%B %d').sort_values().unique()
axs['b'].set_xticklabels(labels=x_dates, rotation=45, ha='right')



# plt.setp(axs['b'].get_xticklabels(), rotation=30, ha='right')
# fig.autofmt_xdate()

sns.despine(trim=True) #offset=10, trim=True)


for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.5, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))


# leaf_gsw.boxplot(column='gsw', ax=axs['b'])
plt.savefig('../figures/leaf_parms.png', dpi=400)


#%%
