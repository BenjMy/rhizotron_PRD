#!/usr/bin/env python
# coding: utf-8

# In[1]:
from pyPRD import processing as proc
import numpy as np
import pyvista as pv
import argparse
import os 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
from scipy.optimize import curve_fit

import datetime
import surveyPRD
from scipy import stats

import locale
import pandas as pd
# Set the locale to English
locale.setlocale(locale.LC_TIME, 'en_US.UTF-8')

#%%

def get_cmd():
    parse = argparse.ArgumentParser()

    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[None], required=False) #[0,1,2,3,4,5]
    period.add_argument('-cycle', '--cycle', nargs='+',
                        # help='list of cycle', type=int, default=[0,1,2,3,4,5,6,7,8,9], required=False)
                        # help='list of cycle', type=int, default=[0,1,2,3,4,5,6,7,8], required=False)
                        help='list of cycle', type=int, default=[0,1,2,3,4,5,6,7,8,9], required=False)
                        # help='list of cycle', type=int, default=[1,2,3,4,5,6,7,8], required=False)
    period.add_argument(
        '-startD', type=str, help='start date to analyse', default=None, required=False)
    period.add_argument('-endD', type=str, help='end date',
                        default=None, required=False)   
    process_param = parse.add_argument_group('process_param')
    process_param.add_argument(
        '-filter_seq', type=int, help='Filter sequence', default=1, required=False)
    process_param.add_argument('-filter_seq_rec', type=int,
                               help='Filter sequence rec', default=0, required=False)
    args = parse.parse_args()
    return(args)

args = get_cmd()

#%% Set path
# -----------------------------------------------------------------------------
imaging_path = '../imaging_ERT_MALM_PaperREV1/'
fmt = '%d/%m/%Y,%H:%M'
startD = None
endD = None


#%% Load log files 
# ----------------------------------------------------------------------------
inversionPathMALM = surveyPRD.definePaths(args)


ERT_log = surveyPRD.load_ERT_survey_log(startDate=args.startD,endDate=args.endD, 
                                        cycles=args.cycle)

# Irrigation schedule
# -------------------
irr_log = surveyPRD.load_irr_log(startDate=args.startD,
                                 endDate=args.endD,
                                 cycles=args.cycle
                                 )

selected_files_ERT = ERT_log[ERT_log['method']=='ERT']['Name'].to_list()
selected_files_MALM = ERT_log[ERT_log['method']=='MALM']['Name'].to_list()


#%% Plot timeline selected window
# ----------------------------------------------------------------------------
ax = surveyPRD.plot_timeline(ERT_log,irr_log,
                        dropMALM=True)
# ax.set_title('Timeline irr. pattern')
plt.savefig('../figures/irr_pattern.png', dpi=400,
            bbox_inches='tight', pad_inches=0)

#%% Load ERT csv data
# -----------------------------------------------------------------------------
file2read=imaging_path+'inversionERT/dfERT' + str([*args.cycle])+'.csv'
df_ERT= pd.read_csv(file2read,index_col=0)

# this is just to get the mesh (kmesh)
kmesh = surveyPRD.process_ERT(imaging_path,
                              selected_files_ERT,
                              ERT_log,
                              reprocessed=0,
                              recip=5,
                              )

#%% Load SWC (Archie) csv data
# -----------------------------------------------------------------------------

file2read=imaging_path+'inversionERT/dfSWC' + str([*args.cycle])+'.csv'
df_SWC= pd.read_csv(file2read,index_col=0)


#%% Load ICSD data
# -----------------------------------------------------------------------------

df_MALM_icsd = pd.read_csv(imaging_path + 
                            inversionPathMALM + 
                            'df_MALM_icsd' + 
                            str([*args.cycle]) 
                            + '.csv',
                            header=[0, 1, 2],
                            parse_dates=True,
                            # index_col=0
                            )

df_MALM_icsd.T.index.names = ('filename','datetime','injType')
df_MALM_icsd = df_MALM_icsd.T

idx = df_MALM_icsd.index
df_MALM_icsd.index = df_MALM_icsd.index.set_levels([idx.levels[0], 
                                                    pd.to_datetime(df_MALM_icsd.index.levels[1]),
                                                    idx.levels[2]
                                                    ]
                                                   )
columns2use = df_MALM_icsd.columns


#%% Load scale data 
# ----------------------------------------------------------------------------
scalesurvey = surveyPRD.load_scale_data(args.startD,args.endD)

scalesurvey= scalesurvey[scalesurvey['datetime']>irr_log['datetime'].min()]


fig, ax = plt.subplots(figsize=(10, 3),constrained_layout=True)
surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax)
ax.set_title('Scale raw data')
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
plt.savefig('../load_cells/figures/scale_raw_data.png', dpi=400)

#%% Pre-process scale data
# ----------------------------------------------------------------------------

scalesurvey=scalesurvey.set_index('datetime')
scalesurvey['cumweight'] = scalesurvey['weight (kg)'].cumsum()
scalesurvey['rolling'] = scalesurvey['weight (kg)'].rolling(12).mean()

scalesurvey['rolling_std'] = scalesurvey['weight (kg)'].rolling(12).std()
scalesurvey['diff'] = scalesurvey['rolling'].diff()

day=[]
scalesurvey.reset_index(inplace=True)
for d in scalesurvey['datetime']:
    if int(d.strftime('%H'))>=7 and int(d.strftime('%H'))<=19:
        day.append(True)
    else:
        day.append(False)

scalesurvey['day'] = day  
scalesurvey['diff']

#%% Plot scale data
# ----------------------------------------------------------------------------

mosaic = '''
            aaaabb
         '''
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=False
                              )

scalesurvey.plot(x='datetime',y='rolling', 
                 ax=axs['a'], 
                 color='k',
                 linestyle='--',
                 label='rolling mean')

scalesurvey.plot.scatter(x='datetime',y='weight (kg)', 
                         color='k',
                         marker='.',
                         ax=axs['a'])
axs['a'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))

axs['a'].spines.right.set_visible(False)
axs['a'].spines.top.set_visible(False)

axin = axs['a'].inset_axes([0.55, 0.72, 0.53, 0.43])

# date0 = '2022-05-20 12:36:00'
# date1 =  '2022-05-25 12:36:00'

date0 = '2022-06-17 12:36:00'
date1 =  '2022-06-20 12:36:00'

scalesurvey_sub = scalesurvey[(scalesurvey['datetime']>date0) & 
                              (scalesurvey['datetime']<date1) 
                              ]
scalesurvey_sub.plot.scatter(x='datetime',y='weight (kg)', 
                         color='k',
                         marker='.',
                         ax=axin)
scalesurvey_sub.plot(x='datetime',y='rolling', 
                 ax=axin, 
                 color='k',
                 label='',
                 linestyle='--',
                 )


axin.fill_between(scalesurvey_sub['datetime'], scalesurvey_sub['rolling'].min(), 
                scalesurvey_sub['rolling'].max(), 
                where=scalesurvey_sub['day']==True,
                color='green', alpha=0.2, label='day'
                )
axin.fill_between(scalesurvey_sub['datetime'], scalesurvey_sub['rolling'].min(),
                scalesurvey_sub['rolling'].max(),
                where=scalesurvey_sub['day']==False,
                color='red', 
                alpha=0.1, 
                label='night'
                )
axin.set_ylabel('') 
axin.set_xlabel('time (Hours)') 
axin.xaxis.set_major_formatter(mdates.DateFormatter('%H'))

axin.xaxis.set_label_position('top') 
axin.xaxis.tick_top()


axs['a'].indicate_inset_zoom(axin)
plt.setp(axs['a'].get_xticklabels(), rotation=30, ha='right')
dt = scalesurvey['datetime'].diff().dt.total_seconds()

scalesurvey['weight_diff']= (scalesurvey['weight (kg)'].diff())/dt
scalesurvey['weight_diff']
df_pos_diff = scalesurvey[scalesurvey['weight_diff'] >=0] # in sec

mean_day= df_pos_diff[df_pos_diff.loc[:,'day']==True]['weight_diff'].mean()
mean_night= df_pos_diff[df_pos_diff.loc[:,'day']==False]['weight_diff'].mean()

index = ['mean']


sec_days = 1.15741e5
rhizo_surface_m2 = 0.03*0.5
kg_m3 = 1e-3
m_mm = 1e-3


mean_day_ax1 = (mean_day*sec_days)
mean_day_ax2 = (mean_day*sec_days*m_mm)*kg_m3/rhizo_surface_m2

mean_night_ax1 = (mean_night*sec_days)
mean_night_ax2 = (mean_night*sec_days*m_mm)*kg_m3/rhizo_surface_m2

df = pd.DataFrame({'day':mean_day_ax1,
                   'night':mean_night_ax1,
                   }, 
                  index=index)
df.plot.bar(rot=0,ylabel=r'$\frac{d_{weight}}{d_{t}}$ (kg/day)',
            ax=axs['b'],
            color=['green','red'],
            alpha=0.1, )

def kgday2mmday(x):
    return (x * m_mm *kg_m3)/(rhizo_surface_m2*1e-6)
    # return x * 2

def mmday2kgday(x):
#     # return x * m_mm *kg_m3/rhizo_surface_m2
#     return x * 3
    return (x * rhizo_surface_m2*1e6)/(m_mm *kg_m3)


secax_y = axs['b'].secondary_yaxis(
    'right', functions=(kgday2mmday, mmday2kgday))
secax_y.set_ylabel(r'(mm/day)')


for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
    
plt.savefig('../figures/day_night_diff.png', dpi=400,bbox_inches='tight', pad_inches=0)


#%% Convert scale variations to soil water content variations
# Finally, from the changes in weight, a change in water content can be calculated. 
# porosity = 0.55
# Rhizotron volume: 7385 cm3
# mixture weight during PRD experiment: 3430 g
# soil_density: 0,46 g/cm3
# density_water = 997 kg/m³

soil_density= (0.46*1e-3)/(1e-6)
density_water = 1000
porosity = 0.55
scalesurvey['weight_diff']

#%% plot PRD effect: ER (Conductivity) VS datetime
# -----------------------------------------------------------------------------

mosaic = '''
            aaaac
            bbbbc
            bbbbd
            bbbbd
         '''
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=False
                              )
irr_log['datetime'].min()
surveyPRD.plot_timeline(ERT_log,irr_log,
                        dropMALM=True,ax=axs['a'],
                        legend=False
                        )
axs['a'].set_xticks([])
axs['a'].set_xlabel('')
ax_ERT = surveyPRD.plot_PRD_effect_ER(kmesh,df_ERT,irr_log=irr_log,
                              unit='mS/m',ax=axs['b'])
axs['b'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))

for label, ax in axs.items():
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))

axs['b'].spines.right.set_visible(False)
axs['b'].spines.top.set_visible(False)
axs['b'].legend(ncol=2, loc='upper center')
axin = axs['c']

date0 = '2022-06-22 00:00:00'
date1 =  '2022-06-24 12:36:00'
axs['b'].add_patch(Rectangle((19165, 0), 4, 60,
                             alpha=0.2,
                             color='lightblue'
                             ),
                   )
axs['c'].add_patch(Rectangle((19165, 25), 4, 35,
                             alpha=0.2,
                             color='lightblue'
                             ),
                   )
ax_ERT = surveyPRD.plot_PRD_effect_ER(kmesh,df_ERT,irr_log=irr_log,
                              unit='mS/m',ax=axin)
ax_ERT.get_legend().remove()
ax_ERT.set_xlim([date0,date1])
axin.xaxis.set_major_locator(mdates.DayLocator(interval=1))
axin.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
axin.set_xlabel('')
axin.set_ylabel('(mS/m)')

# # cycle 6
axin = axs['d']
date0 = '2022-07-05 00:00:00'
date1 =  '2022-07-11 12:36:00'

axs['b'].add_patch(Rectangle((19178, 0), 6, 60,
                             alpha=0.1,
                             color='grey'
                             ),
                   )
axs['b'].get_xlim()
axs['d'].add_patch(Rectangle((19178, 25), 6, 40,
                             alpha=0.1,
                             color='grey'
                             ),
                   )



ax_ERT = surveyPRD.plot_PRD_effect_ER(kmesh,df_ERT,irr_log=irr_log,
                              unit='mS/m',ax=axin)
ax_ERT.get_legend().remove()
ax_ERT.set_xlim([date0,date1])
axin.xaxis.set_major_locator(mdates.DayLocator(interval=2))
axin.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
axin.set_xlabel('')
axin.set_ylabel('(mS/m)')
axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
axs['c'].grid('on', which='major', axis='x',color='0.95' )
axs['c'].grid('on', which='major', axis='y',color='0.95' )
axs['d'].grid('on', which='major', axis='x',color='0.95' )
axs['d'].grid('on', which='major', axis='y',color='0.95' )
axs['a'].set_ylabel('Input water \n (ml)')

plt.savefig('../figures/Cond_variations_datetime.png', 
            dpi=400,
            bbox_inches='tight', 
            pad_inches=0.1
            )

#%% plot PRD_effect: SWC versus datetime
# -----------------------------------------------------------------------------

mosaic = '''
            aaa
            bbb
            bbb
            bbb
         '''
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=True
                              )
irr_log['datetime'].min()
surveyPRD.plot_timeline(ERT_log,irr_log,
                        dropMALM=True,ax=axs['a'],
                        legend=False
                        )
ax_ERT = surveyPRD.plot_PRD_effect_SWC(kmesh,df_SWC,irr_log=irr_log,
                             ax=axs['b'])

axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
    
for label, ax in axs.items():
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
plt.savefig('../figures/SWC_variations_datetime.png', dpi=400)


#%% plot ICSD
# -----------------------------------------------------------------------------
# test

imin = np.loadtxt(imaging_path + 
                    inversionPathMALM + 
                    'vrte_nodes_in_mesh.txt'
                    )

mosaic = '''
            aaa
            bbb
         '''
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=True,
                              figsize=(6,4)
                              )
ax_irr = surveyPRD.plot_timeline(ERT_log,irr_log,
                                 ax=axs['a'],
                                 dropMALM=True,
                                 legend=False)

axs['a'].set_xlabel('')
axs['a'].set_ylabel('Input water \n (ml)'
                    )
                    # , fontsize=8)


_ = surveyPRD.plot_PRD_effect_icsd_gravity(kmesh, 
                                            imin,
                                            df_MALM_icsd, 
                                            irr_log,
                                            ax=axs['b'],
                                            # hours_interval=10,
                                            soil=False,
                                            # color=['red','blue'],
                                            detrend = 0.01
                                            )

axs['b'].set_ylabel('x center \n of mass CSD \n (cm)')
                    # , fontsize=8)
axs['b'].set_ylim([0.20,0.4])
axs['b'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'),
                                   )
plt.setp(axs['b'].get_xticklabels(), rotation=30, ha='right')

axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
axs['b'].spines.right.set_visible(False)
axs['b'].spines.top.set_visible(False)
axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
axs['b'].set_xlabel('date')


plt.savefig('../figures/icsd_PRD_effect.png', dpi=400,
            bbox_inches='tight', pad_inches=0)

 

#%% Compute new statistics for ICSD 
# -----------------------------------------------------------------------------

df_MALM_icsd['medianICSD'] = df_MALM_icsd[columns2use].median(axis=1)
df_MALM_icsd['meanICSD'] = df_MALM_icsd[columns2use].mean(axis=1)
df_MALM_icsd['minICSD'] = df_MALM_icsd[columns2use].min(axis=1)
df_MALM_icsd['maxICSD'] = df_MALM_icsd[columns2use].max(axis=1)
df_MALM_icsd['std'] = df_MALM_icsd[columns2use].std(axis=1)
Q1 = df_MALM_icsd[columns2use].quantile(0.25,axis=1)
Q3 = df_MALM_icsd[columns2use].quantile(0.75,axis=1)
IQR = Q3 - Q1
df_MALM_icsd['IQR'] = IQR
matching8r = [s for i,s in enumerate(list(df_MALM_icsd.index.levels[0])) if "8return" in s]
df_MALM_icsd = df_MALM_icsd.drop(matching8r)



#%% plot petro Archie
# -----------------------------------------------------------------------------

# surveyPRD.petro_plot(
#                         scalesurvey,
#                         kmesh,
#                         df_SWC,
#                         irr_log
#                       )

#%% interpolate weight times on ERT and ICSD times
# -----------------------------------------------------------------------------
# del merge 
Zones=df_ERT['Zone']
df_ERTz = df_ERT.groupby('Zone').mean().T
df_SWCz = df_SWC.groupby('Zone').mean().T
df_SWCz['meanTheta'] = df_SWCz.mean(axis=1)
df_ERTz['meanRes'] = df_ERTz.mean(axis=1)
df_ERTz['meanTheta'] = df_SWCz['meanTheta'] 
df_ERTz.index.name = 'datetime'
df_ERTz.index = pd.to_datetime(df_ERTz.index)
df_ERTz.index.max()
df_ERTz = df_ERTz.reset_index()
scalesurvey.index.min()
scalesurvey.index.max()

#%%
resampled_scalesurvey = scalesurvey.set_index('datetime').resample(rule='1T', 
                                                                   closed='right').nearest()

fig, ax = plt.subplots()
resampled_scalesurvey['weight (kg)'].plot(ax=ax, color='r')
scalesurvey.set_index('datetime')['weight (kg)'].plot(ax=ax)

df_ERTz = df_ERTz.set_index('datetime')
mergeObs=pd.merge(df_ERTz,
                  resampled_scalesurvey,
                  how='left', 
                  on = 'datetime'
                )
# merge['weight (kg)'] = merge['weight (kg)'].interpolate(method ='linear')
mergeObs = mergeObs.reset_index()

# Need to delete time where no scale data were available
# --------------------------------------------------------
mergeObs['weight (kg)'] = mergeObs['weight (kg)'].interpolate(method ='linear')

bool_ind_2del= mergeObs['datetime'].isin([
                                        '2022-05-13 16:25:00',
                                        '2022-07-11 15:50:00',
                                        '2022-07-12 12:00:00'
                        ]
                       )
                        
# mergeObs = mergeObs[~bool_ind_2del]

#%% find nearest date

scalesurvey = scalesurvey.set_index('datetime')
merged_df = pd.merge_asof(df_ERTz, scalesurvey, left_index=True, right_index=True)


#%%
fig, ax = plt.subplots()

ax2 = ax.twinx()  # Create a secondary y-axis

# merged_df['weight (kg)'].plot(ax=ax2, color='k',marker='.')  # Plot the first series on the primary y-axis
mergeObs.set_index('datetime')['weight (kg)'].plot(ax=ax2, color='red',marker='v')  # Plot the second series on the secondary y-axis
df_ERTz['meanRes'].plot(ax=ax, color='blue',marker='.')  # Plot the first series on the primary y-axis
scalesurvey['weight (kg)'].plot(ax=ax2, color='red',marker='v')  # Plot the second series on the secondary y-axis

ax3 = ax.twinx()  # Create a secondary y-axis

df_ERTz['meanTheta'].plot(ax=ax3, color='green',marker='.')  # Plot the first series on the primary y-axis


# Calculate the rolling mean over 4 points for each axis
# rolling_mean_mergeObs = mergeObs.set_index('datetime')['weight (kg)'].rolling(window=6, min_periods=1).mean()
# rolling_mean_df_ERTz = df_ERTz['meanRes'].rolling(window=6, min_periods=1).mean()

# Plot the rolling mean on the corresponding axes
# rolling_mean_mergeObs.plot(ax=ax2, color='red', marker='v', linestyle='--')  # Plot the rolling mean for mergeObs on the secondary y-axis
# rolling_mean_df_ERTz.plot(ax=ax, color='blue', marker='.', linestyle='--')  # Plot the rolling mean for df_ERTz on the primary y-axis

ax2.set_ylabel('Weight (kg)',color='red')  # Label for the primary y-axis
ax.set_ylabel('Mean Res',color='blue')  # Label for the secondary y-axis

plt.show()

#%%
fig, ax = plt.subplots()
ax.scatter(
           mergeObs['meanRes'],
           mergeObs['weight (kg)'],
           )

#%%

#%%
fig, ax = plt.subplots()
ax.scatter(
           mergeObs['meanRes'],
           mergeObs['meanTheta'],
           )



# You are plotting changes in weight over time versus resistivities or versus water contents. 
# Shouldn’t changes in weight be plotted versus changes in water content or changes in resistivity? 
# Although I am not sure whether the latter makes sense since water content and resistivity are non-linearly related in Archie’s Equation. 
# But, what amazes me the most is that resistivity increases when the weight and hence the water content increases. 
# This cannot be correct. Finally, from the changes in weight, a change in water content can be calculated. 
# These changes in water content calculated from weight changes should be one-to-one related to the changes in water content calculated from the ERT measurements. 
# I propose comparing those.



#%%  prepare derivatives 
# -----------------------------------------------------------------------------

PRD4holes = pd.Timestamp('19/05/2022 17:00')
PRD2holes = pd.Timestamp('2022-06-08 11:50')


dt = mergeObs['datetime'].diff().dt.total_seconds()


weight=mergeObs[['datetime','weight (kg)']]
weight['diff'] = weight['weight (kg)'].diff()
weight['deltaT'] = weight['weight (kg)']/dt
weight['diff_deltaT'] = weight['diff']/dt


meanRes = mergeObs[['datetime','meanRes']]
meanRes['diff'] = meanRes['meanRes'].diff()
meanRes['deltaT'] = meanRes['meanRes']/dt
meanRes['diff_deltaT'] = meanRes['diff']/dt


meanTheta = mergeObs[['datetime','meanTheta']]
meanTheta['diff'] = meanTheta['meanTheta'].diff()
meanTheta['deltaT'] = meanTheta['meanTheta']/dt
meanTheta['diff_deltaT'] = meanTheta['diff']/dt



#%% Add some flag (before/after irrigation and before after H4, H2)

before_PRD4holes = []
before_PRD2holes = []

for date in meanTheta.datetime:
    if date <= PRD4holes:
        before_PRD4holes.append(True)
        before_PRD2holes.append(True)
    elif date <= PRD2holes:
        before_PRD4holes.append(False)
        before_PRD2holes.append(True)        
    else:
        before_PRD4holes.append(False)
        before_PRD2holes.append(False)   
     

irr_log

offset_bool_bef_irr = []
offset_bool_after_irr = []
for i, date in enumerate(meanTheta.datetime):
    offset_bef_irr = [date-dirr for dirr in irr_log.datetime]
    offset_bool_bef_irr_tmp = [dirr_delta.total_seconds()>-(60*60*8) and dirr_delta.total_seconds()<0 for dirr_delta in offset_bef_irr]
    offset_bef_irr_insec = [dirr_delta.total_seconds() for dirr_delta in offset_bef_irr]
    if sum(offset_bool_bef_irr_tmp)>0:
        offset_bool_bef_irr.append(True)
    else:
        offset_bool_bef_irr.append(False)


weight['bef_irr_flag'] = offset_bool_bef_irr
weight['before_PRD4holes'] = before_PRD4holes
weight['before_PRD2holes'] = before_PRD2holes


meanRes['bef_irr_flag'] = offset_bool_bef_irr
meanRes['before_PRD4holes'] = before_PRD4holes
meanRes['before_PRD2holes'] = before_PRD2holes

meanTheta['bef_irr_flag'] = offset_bool_bef_irr
meanTheta['before_PRD4holes'] = before_PRD4holes
meanTheta['before_PRD2holes'] = before_PRD2holes


# len(meanTheta)

#%%

# def normalise(p,name):
#     pN=((p[name]-p[name].mean())/p[name].std())
#     return pN

#%%
# coef = meanTheta['diff'].iloc[1]/weight['diff'].iloc[1]
# coef = 1-0.7385
# selec = np.array(before_PRD2holes)

# selec = offset_bool_bef_irr
# selec = ~np.array(offset_bool_bef_irr)
# selec[0] = False

# fig, ax = plt.subplots()
# ax.scatter(
#             weight['diff'][selec]*coef,
#             meanTheta['diff'][selec]
#             )

# ax.scatter(weight['diff'][selec]*coef, meanTheta['diff'][selec])


# z = np.polyfit((weight['diff'][selec]*coef), meanTheta['diff'][selec], 1)
# p = np.poly1d(z)

# _ = ax.plot(weight['diff'][selec]*coef, 
#             p(weight['diff'][selec]*coef), '--')

# ax.set_xlim([-0.1,0.1])
# ax.set_ylim([-0.1,0.1])

# corr_matrix = np.corrcoef(meanTheta['diff'][selec],  
#                           p(weight['diff'][selec]*coef))
# corr = corr_matrix[0,1]
# R_sq = corr**2

# # ax.set_ylabel('Mean Res',color='blue')  # Label for the secondary y-axis
# ax.set_xlabel(r'$\Delta \bar{\Theta}_{scale}$   $[\frac{m^{3}}{m^{3}}]$', fontsize=16)
# ax.set_ylabel(r'$\Delta \bar{\Theta}_{ERT}$   $[\frac{m^{3}}{m^{3}}]$', fontsize=16)
# plt.tight_layout()


#%% ICSD VS RES and THETA variations
# -----------------------------------------------------------------------------
kcolors = ['k']*(len(offset_bool_bef_irr)-1)

threshold = 1e-2 #in % of the total number of source i.e. how much sources carry at least 10% of the current density
df_MALM_icsd.columns
meshXYZ = kmesh[0].mesh.df[['X', 'Y', 'Z']]
xyz_vrte = meshXYZ.iloc[imin]
df_MALM_icsd['threshold'] = (np.sum(df_MALM_icsd[columns2use]>threshold,axis=1)*100)/len(columns2use)

dateStartPRD = df_MALM_icsd.index.get_level_values(1)[0]
df_MALM_icsd_dates_valid = df_MALM_icsd.loc[df_MALM_icsd.index.get_level_values(1)>=dateStartPRD]
beforePRD_nb = len(df_MALM_icsd) - len(df_MALM_icsd_dates_valid)

parm = 'threshold' # medianICSD IQR threshold
pAverage_ICSD_stem = df_MALM_icsd.xs('Stem',level=2)[parm].values 
pAverage_ICSD_soil = df_MALM_icsd.xs('Soil',level=2)[parm].values
len(pAverage_ICSD_stem)

date_icsd_stem = df_MALM_icsd.xs('Stem',level=2).index.get_level_values(1).to_list()
date_icsd_soil = df_MALM_icsd.xs('Soil',level=2).index.get_level_values(1).to_list()
len(date_icsd_stem)

date_icsd_stem_str = [ dd.strftime("%Y-%m-%d") for dd in date_icsd_stem]
date_icsd_soil_str = [ dd.strftime("%Y-%m-%d") for dd in date_icsd_soil]
meanRes_datetime_str = [ dd.strftime("%Y-%m-%d") for dd in meanRes['datetime']]

#%% Merging

df_MALM_icsd_stem = df_MALM_icsd.xs('Stem',level=2)[parm]
df_MALM_icsd_stem = df_MALM_icsd_stem.droplevel('filename')
df_MALM_icsd_soil= df_MALM_icsd.xs('Soil',level=2)[parm]
df_MALM_icsd_soil = df_MALM_icsd_soil.droplevel('filename')

weight = weight.set_index('datetime')
df_MALM_icsd_stem_scale = pd.merge_asof(df_MALM_icsd_stem, weight, left_index=True, right_index=True)

# df_MALM_icsd_soil.sort_index(inplace=True)
df_MALM_icsd_soil_scale = pd.merge_asof(df_MALM_icsd_soil, weight, left_index=True, right_index=True)
df_MALM_icsd_stem_scale = df_MALM_icsd_stem_scale.dropna()
df_MALM_icsd_soil_scale = df_MALM_icsd_soil_scale.dropna()

df_MALM_icsd_soil.index
weight.index

df_MALM_icsd_stem_scale = df_MALM_icsd_stem_scale.loc[df_MALM_icsd_stem_scale.index>=min(df_MALM_icsd_soil_scale.index)]


#%%
# selec = np.array(before_PRD2holes)

mosaic = 'ab'
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharey=False
                              )

for ab in ['a','b']:
    axs[ab].grid('on', which='major', axis='x',color='0.95' )
    axs[ab].grid('on', which='major', axis='y',color='0.95' )
    axs[ab].grid('on', which='major', axis='x',color='0.95' )
    axs[ab].grid('on', which='major', axis='y',color='0.95' )


mergeObs['weight_diff'] = mergeObs['weight (kg)'].diff()
mergeObs['meanTheta_diff'] = mergeObs['meanTheta'].diff()

mergeObs_nonan = mergeObs.dropna()

# coeff = mean(mergeObs['meanTheta_diff']/mergeObs['weight_diff'])

coeff = np.mean((mergeObs['meanTheta_diff']/mergeObs['weight_diff'])[np.isfinite((mergeObs['meanTheta_diff']/mergeObs['weight_diff']))])
# coeff=1
axs['a'].scatter(
            mergeObs_nonan['weight_diff']*coeff,
            mergeObs_nonan['meanTheta_diff'],
            color='k',
            )

# axs['a'].scatter(
#             mergeObs['weight (kg)'][:]*coef,
#             mergeObs['meanTheta'][:],
#             color='k',
#             )

# axs['a'].scatter(weight['diff'][:]*coef, meanTheta['diff'][:])


# z = np.polyfit((mergeObs['weight_diff']*coeff), mergeObs['meanTheta_diff'], 1)
# p = np.poly1d(z)

# _ = ax.plot(mergeObs['weight_diff']*coeff, 
#             p(mergeObs['weight_diff']*coeff), '--')

axs['a'].set_xlim([-0.1,0.1])
axs['a'].set_ylim([-0.1,0.1])

# corr_matrix = np.corrcoef(meanTheta['diff'][:],  
#                           p(weight['diff'][:]*coef))
# corr = corr_matrix[0,1]
# R_sq = corr**2

# mergeObs = mergeObs.dropna()

def func(x, a):
    return a*x


popt, pcov = curve_fit(func, mergeObs_nonan['weight_diff']*coeff, mergeObs_nonan['meanTheta_diff']) # your data x, y to fit
yfit_res_stem = popt[0]*mergeObs_nonan['weight_diff']*coeff
slope, intercept, r_value, p_value, std_err = stats.linregress(mergeObs_nonan['weight_diff']*coeff,
                                                               mergeObs_nonan['meanTheta_diff'],
                                                               )
R2_dres_soil_text = '$R^{2}$=' + str(np.round(r_value,2))
p_dres_soil_text = '$p$=' + str(np.round(p_value,8))


# _ = axs['a'].plot(mergeObs['weight_diff']*coeff, 
#                   yfit_res_stem, '--')


# ax.set_ylabel('Mean Res',color='blue')  # Label for the secondary y-axis
axs['a'].set_xlabel(r'$\Delta \bar{\Theta}_{scale}$   $[\frac{m^{3}}{m^{3}}]$', fontsize=16)
axs['a'].set_ylabel(r'$\Delta \bar{\Theta}_{ERT}$   $[\frac{m^{3}}{m^{3}}]$', fontsize=16)

axs['a'].text(-0.075,0.050,R2_dres_soil_text,c='grey')
axs['a'].text(-0.075,0.0400,p_dres_soil_text,c='grey')
axs['a'].set_aspect('equal')

# Calculate the rolling mean
rolling_mean_stem = df_MALM_icsd_stem_scale['threshold'].rolling(window=8, min_periods=1).mean()
rolling_mean_soil = df_MALM_icsd_soil_scale['threshold'].rolling(window=8, min_periods=1).mean()


df_MALM_icsd_stem_scale.columns

for i, mm in enumerate(['v','+']):
    
    if i==0:
        cond = df_MALM_icsd_stem_scale['bef_irr_flag']==True
        label= 'Ns stem before irr.'
    else:
        cond = df_MALM_icsd_stem_scale['bef_irr_flag']==False
        label= 'Ns stem after irr.'

        
    axs['b'].scatter(x=df_MALM_icsd_stem_scale.index[cond],
                     y=df_MALM_icsd_stem_scale['threshold'].loc[cond],
                     label=label,marker=mm,color='b')
    
axs['b'].scatter(x=df_MALM_icsd_soil_scale.index,
                 y=df_MALM_icsd_soil_scale['threshold'],
                 label='Ns soil',marker='.',color='r')



# Plot the rolling mean
axs['b'].plot(df_MALM_icsd_stem_scale.index, rolling_mean_stem, color='blue', label='Rolling Mean (Stem)',
              linestyle='--')
axs['b'].plot(df_MALM_icsd_soil_scale.index, rolling_mean_soil, color='red', label='Rolling Mean (Soil)',
              linestyle='--')

axs['b'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
plt.setp(axs['b'].get_xticklabels(), rotation=30, ha='right')
axs['b'].set_ylabel('Ns>1% of total CSD (in %)')

axs['b'].legend()
for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

axs['b'].fill_between([min(df_MALM_icsd_soil_scale.index),max(df_MALM_icsd_soil_scale.index)], 2, 10, color='grey', alpha=0.1)
# axs['b'].set_xlim([0.2,0.375])
plt.legend(bbox_to_anchor=(.5, 1.2), loc='center')#, ncol=2)

plt.savefig('../figures/fig9_paper.png', dpi=400)


# axs['a'].scatter(
#                  meanTheta['meanTheta'][bool_date_meanRes_stem],
#                  pAverage_ICSD_stem[bool_date_icsd_stem],
#                  # color='k',
#                  color=kcolors,
#                  marker='v',
#             )



#%%
# #%%
# bool_date_icsd_stem = []
# for date_icsd_i in date_icsd_stem_str:
#     if date_icsd_i in meanRes_datetime_str:
#         bool_date_icsd_stem.append(True)
#     else:
#         bool_date_icsd_stem.append(False)
    

# bool_date_icsd_soil = []
# for date_icsd_i in date_icsd_soil_str:
#     if date_icsd_i in meanRes_datetime_str:
#         bool_date_icsd_soil.append(True)
#     else:
#         bool_date_icsd_soil.append(False)

# bool_date_meanRes_stem = []
# for date_icsd_i in meanRes_datetime_str:
#     if date_icsd_i in date_icsd_stem_str:
#         bool_date_meanRes_stem.append(True)
#     else:
#         bool_date_meanRes_stem.append(False)
        

# bool_date_meanRes_soil = []
# for date_icsd_i in meanRes_datetime_str:
#     if date_icsd_i in date_icsd_soil_str:
#         bool_date_meanRes_soil.append(True)
#     else:
#         bool_date_meanRes_soil.append(False)
# len(bool_date_meanRes_stem)
# len(bool_date_meanRes_soil)

# len(bool_date_icsd_stem)
# len(bool_date_icsd_soil)
# len(meanTheta['meanTheta'])

# #%%
# # in_out_date_meanRes_soil = offset_bool_bef_irr
# # len(offset_bool_bef_irr)
# # len(in_out_date_icsd_soil)

# mosaic = 'ab'
# fig, axs = plt.subplot_mosaic(mosaic,
#                               constrained_layout=True,
#                               sharey=True
#                               )

# axs['a'].grid('on', which='major', axis='x',color='0.95' )
# axs['a'].grid('on', which='major', axis='y',color='0.95' )
# axs['b'].grid('on', which='major', axis='x',color='0.95' )
# axs['b'].grid('on', which='major', axis='y',color='0.95' )


#     # Stem Data
#     # -----------------------------------------------------------------------------
# # len(pAverage_ICSD_stem[in_out_date_icsd_stem])
# # len(meanTheta['meanTheta'])

# axs['a'].scatter(
#                  meanTheta['meanTheta'][bool_date_meanRes_stem],
#                  pAverage_ICSD_stem[bool_date_icsd_stem],
#                  # color='k',
#                  color=kcolors,
#                  marker='v',
#             )


# slope, intercept, r_value, p_value, std_err = stats.linregress(meanTheta['meanTheta'][bool_date_meanRes_stem],
#                                                                pAverage_ICSD_stem[bool_date_icsd_stem])

# # len(meanTheta['meanTheta'][in_out_date_meanRes_soil])
# # len(pAverage_ICSD_soil[in_out_date_icsd_soil])

# axs['a'].scatter(
#                  meanTheta['meanTheta'][bool_date_meanRes_soil],
#                  pAverage_ICSD_soil[bool_date_icsd_soil],
#                  color='blue',
#                  marker='o',
#             )

# yfit_theta_stem = slope*meanTheta['meanTheta'] + intercept
# R2_theta_stem_text = '$R^{2}$=' + str(np.round(r_value,2))
# p_theta_stem_text = '$p$=' + str(np.round(p_value,2))

# # len(weight['weight (kg)'][1:])
# # len(pAverage_ICSD_stem)
# # len(pAverage_ICSD_stem[in_out_date_icsd_stem][1:])
# axs['b'].scatter(
#                   weight['weight (kg)'][1:],
#                   pAverage_ICSD_stem[in_out_date_icsd_stem][1:],
#                   # color='k',
#                   color=kcolors[1:],
#                   label='stem',
#                   marker='v',
#             )

# axs['b'].scatter(
#                   weight['weight (kg)'][in_out_date_meanRes_soil],
#                   pAverage_ICSD_soil[in_out_date_icsd_soil],
#                   color='blue',
#                   label='soil',
#                   marker='o',
#             )



#     # Soil Data
#     # ---------
# len(meanTheta['meanTheta'][in_out_date_meanRes_soil])
# len(weight['diff'])
# len(pAverage_ICSD_soil)

# axs['a'].set_xlabel(r'$\bar{\Theta}$ predicted ($m^{3}/m^{3}$)',fontsize=11)
# axs['b'].set_xlabel(r'$\Delta$ weight (Ohm.m)',fontsize=11) # ($\Omega.m.h^{-1}$')




# # axs['a'].set_ylim([-0.5e-3,1e-3])
# # axs['a'].set_ylabel(parm + ' CSD ($A.m^{-2}$)')
# axs['a'].set_ylabel('Ns>1% of total CSD (in %)')

# for label, ax in axs.items():
#     # label physical distance in and down:
#     trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
#     ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
#             fontsize='medium', verticalalignment='top', fontfamily='serif',
#             bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

# axs['a'].fill_between([0.2,0.375], 2, 10, color='grey', alpha=0.1)
# axs['a'].set_xlim([0.2,0.375])
# plt.legend(bbox_to_anchor=(.5, 1.2), loc='center')
# plt.savefig('../figures/CSD_f_weigth_Res.png', dpi=400)



#%%

if __name__ == '__main__':

    '''
    - Read files
    - Invert ERT
    - Invert MALM
    '''

    # Read and parse args
    # -----------------------------------
    args = get_cmd()
    
    

    