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

import datetime
import surveyPRD
from scipy import stats

#%%

def get_cmd():
    parse = argparse.ArgumentParser()

    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[None], required=False) #[0,1,2,3,4,5]
    period.add_argument('-cycle', '--cycle', nargs='+',
                        # help='list of cycle', type=int, default=[3,4,5,6,7,8], required=False)
                        # help='list of cycle', type=int, default=[0,1,2,3,4,5,6,7,8], required=False)
                        help='list of cycle', type=int, default=[1,2,3,4,5,6,7,8], required=False)
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

#%%
imaging_path = '../imaging_ERT_MALM/'
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
irr_log = surveyPRD.load_irr_log(startDate=args.startD,endDate=args.endD,
                                        cycles=args.cycle)

selected_files_ERT = ERT_log[ERT_log['method']=='ERT']['Name'].to_list()
selected_files_MALM = ERT_log[ERT_log['method']=='MALM']['Name'].to_list()


#%% Plot timeline selected window
# ----------------------------------------------------------------------------
ax = surveyPRD.plot_timeline(ERT_log,irr_log,
                        dropMALM=True)
ax.set_title('Timeline irr. pattern')
plt.savefig('../figures/irr_pattern.png', dpi=400)

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

#%% Load scale data + plot
# ----------------------------------------------------------------------------
scalesurvey = surveyPRD.load_scale_data(args.startD,args.endD)

scalesurvey= scalesurvey[scalesurvey['datetime']>irr_log['datetime'].min()]


fig, ax = plt.subplots(figsize=(10, 3),constrained_layout=True)
surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax)
ax.set_title('Scale raw data')
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
plt.savefig('../load_cells/figures/scale_raw_data.png', dpi=400)

#%%

scalesurvey=scalesurvey.set_index('datetime')
scalesurvey['cumweight'] = scalesurvey['weight (kg)'].cumsum()
#data_scale['rolling'] = data_scale['weight (kg)'].rolling(20).sum()
scalesurvey['rolling'] = scalesurvey['weight (kg)'].rolling(12).mean()

scalesurvey['rolling_std'] = scalesurvey['weight (kg)'].rolling(12).std()
scalesurvey['diff'] = scalesurvey['rolling'].diff()

day=[]
scalesurvey.reset_index(inplace=True)
for d in scalesurvey['datetime']:
    #print(d)
    #print(int(d.strftime('%H')))
    if int(d.strftime('%H'))>=7 and int(d.strftime('%H'))<=19:
        day.append(True)
    else:
        day.append(False)

scalesurvey['day'] = day  
scalesurvey['diff']

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
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
# ax.xaxis.set_minor_locator(mdates.HourLocator(interval=4))

# ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
# ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
axs['a'].xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))

axs['a'].spines.right.set_visible(False)
axs['a'].spines.top.set_visible(False)

# axin = ax.inset_axes([0.55, 0.02, 0.43, 0.43])
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
# axin.set_xticks([])
# axin.set_yticks([]) 

axin.set_ylabel('') 
axin.set_xlabel('time (Hours)') 

# axin.xaxis.set_major_locator(mdates.DayLocator(interval=3))
# axin.xaxis.set_minor_locator(mdates.HourLocator(interval=6))
# axin.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
axin.xaxis.set_major_formatter(mdates.DateFormatter('%H'))

axin.xaxis.set_label_position('top') 
axin.xaxis.tick_top()


axs['a'].indicate_inset_zoom(axin)
# axs['a'].autofmt_xdate()
# plt.setp(axs['a'].xaxis.xticks()[1], rotation=30, ha='right') # ha is the same as horizontalalignment

plt.setp(axs['a'].get_xticklabels(), rotation=30, ha='right')
# delete all positive variations
# cumsum during diurn/night

# scalesurvey['weight_diff']= scalesurvey['weight (kg)'].diff()

dt = scalesurvey['datetime'].diff().dt.total_seconds()

scalesurvey['weight_diff']= (scalesurvey['weight (kg)'].diff())/dt
scalesurvey['weight_diff']

df_pos_diff = scalesurvey[scalesurvey['weight_diff'] >=0] # in sec
#df_pos_diff = df[df['weight_diff'] >0]

mean_day= df_pos_diff[df_pos_diff.loc[:,'day']==True]['weight_diff'].mean()
mean_night= df_pos_diff[df_pos_diff.loc[:,'day']==False]['weight_diff'].mean()

index = ['mean']
#*5*1e3/60


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
#ax = df.plot.bar(rot=0,ylabel=mean ')
df.plot.bar(rot=0,ylabel=r'$\frac{d_{weight}}{d_{t}}$ (kg/day)',
            ax=axs['b'],
            color=['green','red'],
            alpha=0.1, )

# def deg2rad(x):
#     return x * np.pi / 180


# def rad2deg(x):
#     return x * 180 / np.pi


# def kgday2mmh(x):
#     return (x * m_mm *kg_m3)/rhizo_surface_m2
#     # return x * 2

# def mmh2kgday(x):
# #     # return x * m_mm *kg_m3/rhizo_surface_m2
# #     return x * 3
#     return (x * rhizo_surface_m2)/(m_mm *kg_m3)

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
# 
#%%




#df['ET_sec'] = (df_pos_diff['weight_diff']*5*1e3/60)*rhizo_surface_m2
#data_scale['ET_h'] =  (weight_lost/round(delta.total_seconds()/60))*rhizo_surface_m2


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

#%%
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

# df_MALM_icsd.columns
# df_MALM_icsd.mode(axis=1).mean(axis=1).groupby('injType').plot()
# df_MALM_icsd.index.levels[1]

# df_MALM_icsd = df_MALM_icsd.T.reset_index()
# df_MALM_icsd['datetime'] = df_MALM_icsd['datetime'].astype('datetime64[ns]')
# df_MALM_icsd.set_index(['filename', 'datetime', 'injection'], inplace=True)
# df_MALM_icsd = df_MALM_icsd.T


imin = np.loadtxt(imaging_path + 
                    inversionPathMALM + 
                    'vrte_nodes_in_mesh.txt'
                    )
# len(imin)

# df_MALM_icsd['medianICSD'].groupby('injType').mean()

#%% plot PRD_effect: ER
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
                        dropMALM=True,ax=axs['a'])

axs['a'].set_xticks([])
axs['a'].set_xlabel('')

ax_ERT = surveyPRD.plot_PRD_effect_ER(kmesh,df_ERT,irr_log=irr_log,
                              unit='mS/m',ax=axs['b'])
for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))

axs['b'].spines.right.set_visible(False)
axs['b'].spines.top.set_visible(False)

# axin = axs['b'].inset_axes([0.55, 0.98, 0.23, 0.53])
# # cycle 4
axin = axs['c']

date0 = '2022-06-22 00:00:00'
date1 =  '2022-06-24 12:36:00'
ax_ERT = surveyPRD.plot_PRD_effect_ER(kmesh,df_ERT,irr_log=irr_log,
                              unit='mS/m',ax=axin)
ax_ERT.get_legend().remove()
ax_ERT.set_xlim([date0,date1])
# axin.set_ylabel('')
# axin.set_xlabel('')
# axin.xaxis.set_label_position('top') 
# axin.xaxis.tick_top()
# plt.setp(axin.get_xticklabels(), rotation=30, ha='right')
axin.xaxis.set_major_locator(mdates.DayLocator(interval=1))
axin.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
# axs['b'].indicate_inset_zoom(axin)

axin.set_xlabel('')
axin.set_ylabel('(mS/m)')


# axin = axs['b'].inset_axes([0.95, 0.98, 0.23, 0.53])
# # cycle 6
axin = axs['d']
date0 = '2022-07-05 00:00:00'
date1 =  '2022-07-11 12:36:00'
ax_ERT = surveyPRD.plot_PRD_effect_ER(kmesh,df_ERT,irr_log=irr_log,
                              unit='mS/m',ax=axin)
ax_ERT.get_legend().remove()
ax_ERT.set_xlim([date0,date1])
# axin.set_ylabel('')
# axin.set_xlabel('')
# axin.xaxis.set_label_position('top') 
# axin.xaxis.tick_top()
# plt.setp(axin.get_xticklabels(), rotation=30, ha='right')
axin.xaxis.set_major_locator(mdates.DayLocator(interval=2))
axin.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
# axs['b'].indicate_inset_zoom(axin)

axin.set_xlabel('')
axin.set_ylabel('(mS/m)')

# from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

# axins = zoomed_inset_axes(axs['b'], 6, loc=1) 

# for axis in ['top','bottom','left','right']:
#     axin.spines[axis].set_linewidth(3)
#     axin.spines[axis].set_color('r')
    
# mark_inset(axs['b'], axin, loc1=2, loc2=1, fc="none", lw=1, ec='k')


axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
axs['c'].grid('on', which='major', axis='x',color='0.95' )
axs['c'].grid('on', which='major', axis='y',color='0.95' )
axs['d'].grid('on', which='major', axis='x',color='0.95' )
axs['d'].grid('on', which='major', axis='y',color='0.95' )
    
axs['a'].set_ylabel('Input water \n (mL)')

# 29/6/2022 (cycle 5)
# Day + 6 (5/7/2022)

plt.savefig('../figures/Cond_variations.png', dpi=400,bbox_inches='tight', pad_inches=0)




#%% plot PRD_effect: SWC
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
                        dropMALM=True,ax=axs['a'])

df_ERT
ax_ERT = surveyPRD.plot_PRD_effect_SWC(kmesh,df_SWC,irr_log=irr_log,
                             ax=axs['b'])

axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
    

for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
# ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
plt.savefig('../figures/Cond_variations.png', dpi=400)


 
#%% plot petro Archie
# -----------------------------------------------------------------------------

# surveyPRD.petro_plot(
#                         scalesurvey,
#                         kmesh,
#                         df_SWC,
#                         irr_log
#                       )

#%% interpolate weight times on ERT times

scalesurvey
Zones=df_ERT['Zone']


df_ERTz = df_ERT.groupby('Zone').mean().T
df_SWCz = df_SWC.groupby('Zone').mean().T
df_SWCz['meanTheta'] = df_SWCz.mean(axis=1)
df_ERTz['meanRes'] = df_ERTz.mean(axis=1)
df_ERTz['meanTheta'] = df_SWCz['meanTheta'] 

df_ERTz.index.name = 'datetime'
df_ERTz.index = pd.to_datetime(df_ERTz.index)

df_ERTz.index.max()

scalesurvey.index.min()
scalesurvey.index.max()

merge=pd.merge(df_ERTz,scalesurvey,
               how='outer', 
                left_on='datetime', 
                right_on='datetime',
                sort=True,
               )
merge.columns
merge['weight (kg)'] = merge['weight (kg)'].interpolate(method ='linear')

merge=merge[pd.isna(merge['Left'])==False]

#%% WEIGHT VS RES and THETA variations

# prepare derivatives 

dt = merge['datetime'].diff().dt.total_seconds()
weightN_diff=merge['weight (kg)']/dt
meanResN_diff = merge['meanRes'].diff()/ dt
meanThetaN_diff = merge['meanTheta'].diff()/dt
# meanThetaN_diff = merge['meanTheta']


mosaic = 'ab'
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=True
                              )

x11 = np.linspace(np.nanmin(weightN_diff),np.nanmax(weightN_diff),10)
axs['a'].scatter(weightN_diff,
            meanResN_diff,
            color='k',
            )
axs['a'].set_xlabel(r'$\frac{d_{weight}}{dt}$', fontsize=14)
axs['a'].set_ylabel(r'$\frac{d \bar{\rho}}{dt}$', fontsize=16)

axs['b'].scatter(weightN_diff,
            meanThetaN_diff,
            color='k',
            )
axs['b'].set_xlabel(r'$\frac{d_{weight}}{dt}$', fontsize=14)
axs['b'].set_ylabel(r'$\frac{d \bar{\Theta}}{dt}$', fontsize=14)

slope, intercept, r_value, p_value, std_err = stats.linregress(weightN_diff[1:], meanThetaN_diff[1:])
yfit_theta = slope*weightN_diff[1:] + intercept
R2_theta_text = '$R^{2}$=' + '{:.3f}'.format(r_value)
p_theta_text = '$p$=' + '{:.1e}'.format(p_value)


slope, intercept, r_value, p_value, std_err = stats.linregress(weightN_diff[1:], meanResN_diff[1:])
yfit_res = slope*weightN_diff[1:] + intercept
R2_res_text = '$R^{2}$=' + '{:1.3f}'.format(r_value)
p_res_text = '$p$=' + '{:.1e}'.format(p_value)

axs['a'].plot(weightN_diff[1:],yfit_res,
          linestyle='--',
          linewidth=0.1,
          color='k')

axs['b'].plot(weightN_diff[1:],yfit_theta,
          linestyle='--',
          linewidth=0.1,
          color='k')

axs['a'].text(0,max(yfit_res)/2,R2_res_text)
axs['a'].text(0,max(yfit_res)/2-2e-4,p_res_text)
axs['b'].text(0,max(yfit_theta)/2,R2_theta_text)
axs['b'].text(0,max(yfit_theta)/2-1e-6,p_theta_text)


for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))


plt.savefig('../figures/weight_f_theta_Res.png', dpi=400)


#%%

def normalise(p,name):
    pN=((p[name]-p[name].mean())/p[name].std())
    return pN

#%%

# img=np.random.randint(0,2,(100,100)).astype(np.uint8)

# x = range(0, img.shape[0])
# y = range(0, img.shape[1])

# (X,Y) = np.meshgrid(x,y)

# x_coord = (X*img).sum() / img.sum().astype("float")
# y_coord = (Y*img).sum() / img.sum().astype("float")

# plt.imshow(img)

#%% ICSD VS RES and THETA variations

# threshold = 1e-2 in % of the total number of source i.e. how much sources carry at least 10% of the current density
df_MALM_icsd.columns

meshXYZ = kmesh[0].mesh.df[['X', 'Y', 'Z']]
xyz_vrte = meshXYZ.iloc[imin]

# df_MALM_icsd[]

# df_MALM_icsd.
# 8/6/2022	12:30
df_MALM_icsd['threshold'] = (np.sum(df_MALM_icsd[columns2use]>1e-2,axis=1)*100)/len(columns2use)

dateStartPRD = '2022-06-08 12:30:00'
# dateStartPRD = '2022-06-05 12:30:00'
df_MALM_icsd_dates_valid = df_MALM_icsd.loc[df_MALM_icsd.index.get_level_values(1)>dateStartPRD]
# pAverage_ICSD_stem = df_MALM_icsd.xs((dates_valid,'Stem'),level=[0,2])[parm].values 
# df_MALM_icsd_dates_valid['datetime']
beforePRD_nb = len(df_MALM_icsd) - len(df_MALM_icsd_dates_valid)

meanRes = df_ERTz['meanRes']
meanTheta = df_ERTz['meanTheta']
dt = (meanRes.reset_index()['datetime'].diff().dt.total_seconds()).values

parm = 'threshold' # medianICSD IQR threshold
pAverage_ICSD_stem = df_MALM_icsd.xs('Stem',level=2)[parm].values 
pAverage_ICSD_soil = df_MALM_icsd.xs('Soil',level=2)[parm].values

# id2keep = list(pAverage_ICSD_stem<5e9)


# df_MALM_icsd['threshold']

dmeanRes_dt = meanRes.diff()[1:]/(meanRes.reset_index()['datetime'].diff().dt.total_seconds()).values[1:]

# id2keep = list(
#                 np.arange(int(sum(df_MALM_icsd.index.get_level_values(0)<dateStartPRD)/2)-1,
#                          len(dmeanRes_dt)
#                          )
#                )

id2keep = list(np.arange(0,len(pAverage_ICSD_soil)))

# len(id2keep)
mosaic = 'ab'
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharey=True
                              )

axs['a'].grid('on', which='major', axis='x',color='0.95' )
axs['a'].grid('on', which='major', axis='y',color='0.95' )
axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )

# Stem Data
# ---------
len(meanTheta[id2keep])
len(meanTheta)
len(pAverage_ICSD_stem)
len(pAverage_ICSD_soil)
len(id2keep)
len(dmeanRes_dt)
axs['a'].scatter(
                 meanTheta[id2keep],
                 pAverage_ICSD_stem,
                 color='k',
                 marker='v',
            )

axs['a'].scatter(
                 meanTheta[0:beforePRD_nb],
                 pAverage_ICSD_stem[0:beforePRD_nb],
                 color='grey',
                 marker='v',
            )


slope, intercept, r_value, p_value, std_err = stats.linregress(meanTheta[id2keep],pAverage_ICSD_stem[id2keep])
yfit_theta_stem = slope*meanTheta + intercept
R2_theta_stem_text = '$R^{2}$=' + str(np.round(r_value,2))
p_theta_stem_text = '$p$=' + str(np.round(p_value,2))


axs['b'].scatter(
                  dmeanRes_dt[:],
                  pAverage_ICSD_stem[1:],
                  color='k',
                  label='stem',
                  marker='v',
            )


axs['b'].scatter(
                  dmeanRes_dt[0:beforePRD_nb],
                  pAverage_ICSD_stem[1:beforePRD_nb+1],
                  color='grey',
                  label='stem (before PRD)',
                  marker='v',
            )


slope, intercept, r_value, p_value, std_err = stats.linregress(dmeanRes_dt[:],pAverage_ICSD_stem[1:])
yfit_meanRes_stem = slope*dmeanRes_dt + intercept
R2_meanRes_stem_text = '$R^{2}$=' + str(np.round(r_value,2))
p_meanRes_stem_text = '$p$=' + str(np.round(p_value,2))

# Soil Data
# ---------

axs['a'].scatter(
                 meanTheta[id2keep],
                 pAverage_ICSD_soil,
                 color='blue',
                 marker='o',
            )

axs['b'].scatter(
                  dmeanRes_dt[:],
                  pAverage_ICSD_soil[1:],
                  color='blue',
                  label='soil',
                  marker='o',
            )


axs['a'].set_xlabel(r'$\Theta$ predicted ($m^{3}/m^{3}$)')
axs['b'].set_xlabel(r'$\frac{d\rho}{dt}$',fontsize=14) # ($\Omega.m.h^{-1}$')




# axs['a'].set_ylim([-0.5e-3,1e-3])
# axs['a'].set_ylabel(parm + ' CSD ($A.m^{-2}$)')
axs['a'].set_ylabel('Ns>1% of total CSD (in %)')

for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

axs['a'].fill_between([0.2,0.375], 2, 10, color='grey', alpha=0.1)
axs['a'].set_xlim([0.2,0.375])
plt.legend(bbox_to_anchor=(.5, 1.2), loc='center')
plt.savefig('../figures/CSD_f_weigth_Res.png', dpi=400)

#%%

#%% ICSD VS RES and THETA variations

# threshold = 1e-2 in % of the total number of source i.e. how much sources carry at least 10% of the current density
df_MALM_icsd.columns

meshXYZ = kmesh[0].mesh.df[['X', 'Y', 'Z']]
xyz_vrte = meshXYZ.iloc[imin]

# df_MALM_icsd[]

# df_MALM_icsd.
# 8/6/2022	12:30
df_MALM_icsd['threshold'] = (np.sum(df_MALM_icsd[columns2use]>1e-2,axis=1)*100)/len(columns2use)

dateStartPRD = '2022-06-08 12:30:00'

dt = (meanRes.reset_index()['datetime'].diff().dt.total_seconds())
# dateStartPRD = '2022-06-05 12:30:00'
df_MALM_icsd_dates_valid = df_MALM_icsd.loc[df_MALM_icsd.index.get_level_values(1)>dateStartPRD]
# pAverage_ICSD_stem = df_MALM_icsd.xs((dates_valid,'Stem'),level=[0,2])[parm].values 
# df_MALM_icsd_dates_valid['datetime']
beforePRD_nb = len(df_MALM_icsd) - len(df_MALM_icsd_dates_valid)

meanRes = df_ERTz['meanRes']
meanTheta = df_ERTz['meanTheta']
dt = (meanRes.reset_index()['datetime'].diff().dt.total_seconds()).values

meanTheta_diff = meanTheta.diff()/dt
parm = 'threshold' # medianICSD IQR threshold
pAverage_ICSD_stem = df_MALM_icsd.xs('Stem',level=2)[parm].values 
pAverage_ICSD_soil = df_MALM_icsd.xs('Soil',level=2)[parm].values

pAverage_ICSD_stem_diff = df_MALM_icsd.xs('Stem',level=2)[parm].diff().values/dt
pAverage_ICSD_soil_diff = df_MALM_icsd.xs('Soil',level=2)[parm].diff().values/dt



# id2keep = list(pAverage_ICSD_stem<5e9)


# df_MALM_icsd['threshold']

dmeanRes_dt = meanRes.diff()[1:]/dt[1:]

# id2keep = list(
#                 np.arange(int(sum(df_MALM_icsd.index.get_level_values(0)<dateStartPRD)/2)-1,
#                          len(dmeanRes_dt)
#                          )
#                )

id2keep = list(np.arange(0,len(pAverage_ICSD_soil)))

# len(id2keep)
mosaic = 'ab'
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharey=True
                              )

# Stem Data
# ---------
axs['a'].scatter(
                 meanTheta_diff,
                 pAverage_ICSD_stem_diff,
                 color='k',
                 marker='v',
            )

axs['a'].scatter(
                 meanTheta_diff[0:beforePRD_nb],
                 pAverage_ICSD_stem_diff[0:beforePRD_nb],
                 color='grey',
                 marker='v',
            )


slope, intercept, r_value, p_value, std_err = stats.linregress(meanTheta[id2keep],pAverage_ICSD_stem[id2keep])
yfit_theta_stem = slope*meanTheta + intercept
R2_theta_stem_text = '$R^{2}$=' + str(np.round(r_value,2))
p_theta_stem_text = '$p$=' + str(np.round(p_value,2))


axs['b'].scatter(
                  dmeanRes_dt[:],
                  pAverage_ICSD_stem_diff[1:],
                  color='k',
                  label='stem',
                  marker='v',
            )


axs['b'].scatter(
                  dmeanRes_dt[0:beforePRD_nb],
                  pAverage_ICSD_stem_diff[1:beforePRD_nb+1],
                  color='grey',
                  label='stem (before PRD)',
                  marker='v',
            )


slope, intercept, r_value, p_value, std_err = stats.linregress(dmeanRes_dt[:],pAverage_ICSD_stem[1:])
yfit_meanRes_stem = slope*dmeanRes_dt + intercept
R2_meanRes_stem_text = '$R^{2}$=' + str(np.round(r_value,2))
p_meanRes_stem_text = '$p$=' + str(np.round(p_value,2))

# Soil Data
# ---------

axs['a'].scatter(
                 meanTheta_diff,
                 pAverage_ICSD_soil_diff,
                 color='blue',
                 marker='o',
            )

axs['b'].scatter(
                  dmeanRes_dt[:],
                  pAverage_ICSD_soil_diff[1:],
                  color='blue',
                  label='soil',
                  marker='o',
            )


# axs['a'].set_xlabel(r'd$\Theta$ predicted ($m^{3}/m^{3}$)')
axs['a'].set_xlabel(r'$\frac{d\Theta}{dt}$ ($m^{3}/m^{3}$)')
axs['b'].set_xlabel(r'$\frac{d\rho}{dt}$',fontsize=14) # ($\Omega.m.h^{-1}$')


axs['a'].grid('on', which='major', axis='x',color='0.95' )
axs['a'].grid('on', which='major', axis='y',color='0.95' )
axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )

# axs['a'].set_ylim([-0.5e-3,1e-3])
# axs['a'].set_ylabel(parm + ' CSD ($A.m^{-2}$)')

# axs['a'].set_ylabel('d(Ns>1% of total CSD (in %))')
axs['a'].set_ylabel(r'$\frac{dNs}{dt}$')

# axs['a'].set_ylabel('d(Ns>1% of total CSD (in %))')
# axs['a'].set_ylabel('Ns>1% of total CSD (in %)')

for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

# axs['a'].fill_between([0.2,0.375], 2, 10, color='grey', alpha=0.1)
# axs['a'].set_xlim([0.2,0.375])
plt.legend(bbox_to_anchor=(.5, 1.2), loc='center')
plt.savefig('../figures/CSD_f_weigth_Res4.png', dpi=400)


#%%
grp = df_MALM_icsd['threshold'].groupby(by='injType')
# grp = grp.droplevel(level=[0,2])

# grp.plot.scatter(x='datetime',y='threshold')

#%% DEPRECATED

# mergeICSD_SWC=pd.merge(df_MALM_icsd,df_ERTz,
#                # how='outer', 
#                right_on=['datetime','injType'],
#                left_index=True,
#                sort=True,
#                )





# # len(df_MALM_icsd)
# len(id2keep)
# len(meanTheta)

# medianICSD_stemN = normalise(df_MALM_icsd.xs('Stem',level=2),'medianICSD')
# medianICSD_soilN = normalise(df_MALM_icsd.xs('Soil',level=2),'medianICSD')

# meanResN = normalise(df_ERTz,'meanRes')

# fig, ax = plt.subplots(1,2)

# x11 = np.arange(-2,3)
# axs['a'].plot(x11,x11,
#          linestyle='--',
#          color='k')
# axs['a'].scatter(
#                  mergeICSD_SWC['meanTheta'],
#                  mergeICSD_SWC['meanRes'],
#                  color='k',
#             )


mosaic = 'ab'
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharey=True
                              )

# len(medianICSD_stem[:])
# len(meanTheta)

axs['a'].scatter(
                 meanTheta[id2keep],
                 medianICSD_stem[id2keep],
                 color='k',
                  marker='v',
            )

popt, pcov = curve_fit(f, meanTheta, medianICSD_stem) # your data x, y to fit
yfit_theta_stem = popt[0]*meanTheta + popt[1]

slope, intercept, r_value, p_value, std_err = stats.linregress(meanTheta[id2keep],medianICSD_stem[id2keep])
yfit_theta_stem = slope*meanTheta + intercept
R2_theta_stem_text = '$R^{2}$=' + str(np.round(r_value,2))
p_theta_stem_text = '$p$=' + str(np.round(p_value,2))

    
axs['a'].plot(
                  meanTheta[id2keep],
                  yfit_theta_stem[id2keep],
                  color='k',
                  linestyle='--',
            )


axs['a'].scatter(
                 meanTheta[id2keep],
                 medianICSD_soil[id2keep],
                 color='grey',
                  marker='o',
            )

popt, pcov = curve_fit(f, meanTheta, medianICSD_soil) # your data x, y to fit
yfit_theta_soil = popt[0]*meanTheta + popt[1]

slope, intercept, r_value, p_value, std_err = stats.linregress(meanTheta,medianICSD_soil)
yfit_theta_soil = slope*meanTheta + intercept
R2_theta_soil_text = '$R^{2}$=' + str(np.round(r_value,2))
p_theta_soil_text = '$p$=' + str(np.round(p_value,2))


axs['a'].text(0.25,2.5e-3,R2_theta_soil_text,c='grey')
axs['a'].text(0.25,2.3e-3,p_theta_soil_text,c='grey')


axs['a'].text(0.30,2.5e-3,R2_theta_stem_text)
axs['a'].text(0.30,2.3e-3,p_theta_stem_text)



axs['a'].plot(
                  meanTheta[id2keep],
                  yfit_theta_soil[id2keep],
                  color='grey',
                  linestyle='--'
            )


axs['b'].scatter(
                  drho_dt[:],
                  medianICSD_stem[1:],
                  color='k',
                  label='stem',
                  marker='v',
            )

popt, pcov = curve_fit(f, meanRes.diff()[1:], medianICSD_stem[1:]) # your data x, y to fit
yfit_res_stem = popt[0]*meanRes.diff()[1:] + popt[1]
slope, intercept, r_value, p_value, std_err = stats.linregress(meanRes.diff()[1:],medianICSD_stem[1:])
yfit_dres_soil = slope*meanRes.diff()[1:] + intercept
R2_dres_soil_text = '$R^{2}$=' + str(np.round(r_value,2))
p_dres_soil_text = '$p$=' + str(np.round(p_value,2))


# axs['b'].plot(
#                   drho_dt[:],
#                   yfit_res_stem,
#                   color='k',
#                   linestyle='-',
#             )

# axs['b'].scatter(
#                   drho_dt[:],
#                   medianICSD_soil[1:],
#                   color='grey',
#                   label='soil',
#                   marker='o',
#                   # s=10,
#             )


# len(medianICSD_soil)
# len(drho_dt)

popt, pcov = curve_fit(f, meanRes.diff()[1:], medianICSD_soil[1:]) # your data x, y to fit
yfit_res = popt[0]*meanRes.diff()[1:] + popt[1]
slope, intercept, r_value, p_value, std_err = stats.linregress(meanRes.diff()[1:],medianICSD_soil[1:])
yfit_dres_stem = slope*meanRes.diff()[1:] + intercept
R2_dres_stem_text = '$R^{2}$=' + str(np.round(r_value,2))
p_dres_stem_text = '$p$=' + str(np.round(p_value,2))


# axs['b'].text(-10,2.5e-3,R2_dres_soil_text,c='grey')
# axs['b'].text(-10,2.3e-3,p_dres_soil_text,c='grey')


# axs['b'].text(0,2.5e-3,R2_dres_stem_text)
# axs['b'].text(0,2.3e-3,p_dres_stem_text)



# axs['b'].plot(
#                   drho_dt[1:],
#                   yfit_res,
#                   color='grey',
#                   linestyle='--',
#                   # linewidth=1,
#                   # marker='+'
#             )


                    # style=['+-', 'o-', '.--', 's:']   


# axs['c'].scatter(
#                   medianICSD_stem,
#                   meanTheta,
#                   color='k',
#             )

# axs['c'].scatter(
#                   medianICSD_soil,
#                   meanTheta,
#                   color='b',
#             )

axs['b'].legend()
axs['a'].set_ylabel('Median CSD')
# axs['b'].set_ylabel('Standart Dev. CSD')



# axs['a'].set_xlabel('Res')
axs['a'].set_xlabel(r'$\Theta$ ($m^{3}/m^{3}$)')
axs['b'].set_xlabel(r'$\frac{d\rho}{dt}$',fontsize=14) # ($\Omega.m.h^{-1}$')
# axs['c'].set_xlabel('$d_{weight}/dt$')

axs['a'].set_ylim([-0.5e-3,1e-3])
axs['a'].set_ylabel('Median CSD ($A.m^{-2}$)')

# axs['a'].axis('square')
# axs['b'].axis('square')
plt.savefig('../figures/CSD_f_weigth_Res.png', dpi=400)

#%% plot ICSD
# -----------------------------------------------------------------------------

mosaic = '''
            aaac
            bbbd
         '''
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=False
                              )


# fig, axs = plt.subplots(2,sharex=True)
# color = 'tab:red'
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y1-axis', color = color)
# ax2 = ax.twinx()
# color = 'tab:green'
# ax2.set_ylabel('Y2-axis', color = color)     
# fig, axs = plt.subplots(3,sharex=False)
ax_irr = surveyPRD.plot_timeline(ERT_log,irr_log,ax=axs['a'],dropMALM=True)
# ax_scale = surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax, color='red')
# ax_SWC = surveyPRD.plot_PRD_effect_SWC(kmesh,df_SWC,irr_log=irr_log,ax=ax2,
#                                        LR=True, topdown=False,
#                                        color=['green','yellow'])

df_MALM_icsdT = df_MALM_icsd.drop('medianICSD',axis=1)
_ = surveyPRD.plot_PRD_effect_icsd(kmesh, 
                                    imin,
                                    df_MALM_icsdT, irr_log,
                                    ax=axs['b'],
                                    hours_interval=10,
                                    soil=False,
                                    # color=['red','blue'],
                                    detrend = 0.05
                                    )
# plt.legend(['SWC Left','SWC Right','ICSD Left','ICSD right'])
# ax_ERT = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,irr_log=irr_log,ax=ax2)
# surveyPRD.petro_plots(scalesurvey,k_indiv_merged,df_SWC,irr_log,ax=ax2)
# ax.set_title('icsd VS weight data')
axs['b'].set_ylabel('CSD')

ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )


axs['b'].spines.right.set_visible(False)
axs['b'].spines.top.set_visible(False)

# axin = axs['b'].inset_axes([0.55, 0.98, 0.23, 0.53])
# # cycle 4
axin = axs['c']

date0 = '2022-06-22 00:00:00'
date1 =  '2022-06-24 12:36:00'
ax_ERT = surveyPRD.plot_PRD_effect_icsd(kmesh, 
                                    imin,
                                    df_MALM_icsdT, irr_log,
                                    ax=axin,
                                    hours_interval=10,
                                    soil=False,
                                    # color=['red','blue'],
                                    detrend = 0.05
                                    )
ax_ERT.get_legend().remove()
ax_ERT.set_xlim([date0,date1])
# axin.set_ylabel('')
# axin.set_xlabel('')
# axin.xaxis.set_label_position('top') 
# axin.xaxis.tick_top()
# plt.setp(axin.get_xticklabels(), rotation=30, ha='right')
axin.xaxis.set_major_locator(mdates.DayLocator(interval=1))
axin.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
# axs['b'].indicate_inset_zoom(axin)

axin.set_xlabel('')
axin.set_ylabel('CSD')


# axin = axs['b'].inset_axes([0.95, 0.98, 0.23, 0.53])
# # cycle 6
axin = axs['d']
date0 = '2022-07-05 00:00:00'
date1 =  '2022-07-11 12:36:00'
ax_ERT = surveyPRD.plot_PRD_effect_icsd(kmesh, 
                                    imin,
                                    df_MALM_icsdT, irr_log,
                                    ax=axin,
                                    hours_interval=10,
                                    soil=False,
                                    # color=['red','blue'],
                                    detrend = 0.05
                                    )
ax_ERT.get_legend().remove()
ax_ERT.set_xlim([date0,date1])
# axin.set_ylabel('')
# axin.set_xlabel('')
# axin.xaxis.set_label_position('top') 
# axin.xaxis.tick_top()
# plt.setp(axin.get_xticklabels(), rotation=30, ha='right')
axin.xaxis.set_major_locator(mdates.DayLocator(interval=2))
axin.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
# axs['b'].indicate_inset_zoom(axin)

axin.set_xlabel('')
axin.set_ylabel('CSD')

# from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

# axins = zoomed_inset_axes(axs['b'], 6, loc=1) 

# for axis in ['top','bottom','left','right']:
#     axin.spines[axis].set_linewidth(3)
#     axin.spines[axis].set_color('r')
    
# mark_inset(axs['b'], axin, loc1=2, loc2=1, fc="none", lw=1, ec='k')


axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
axs['c'].grid('on', which='major', axis='x',color='0.95' )
axs['c'].grid('on', which='major', axis='y',color='0.95' )
axs['d'].grid('on', which='major', axis='x',color='0.95' )
axs['d'].grid('on', which='major', axis='y',color='0.95' )


plt.savefig('../figures/icsd_PRD_effect.png', dpi=400)

#%% plot ICSD
# -----------------------------------------------------------------------------

mosaic = '''
            aaac
            bbbd
         '''
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=True
                              )


# fig, axs = plt.subplots(2,sharex=True)
# color = 'tab:red'
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y1-axis', color = color)
# ax2 = ax.twinx()
# color = 'tab:green'
# ax2.set_ylabel('Y2-axis', color = color)     
# fig, axs = plt.subplots(3,sharex=False)
ax_irr = surveyPRD.plot_timeline(ERT_log,irr_log,ax=axs['a'],dropMALM=True)
# ax_scale = surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax, color='red')
# ax_SWC = surveyPRD.plot_PRD_effect_SWC(kmesh,df_SWC,irr_log=irr_log,ax=ax2,
#                                        LR=True, topdown=False,
#                                        color=['green','yellow'])

df_MALM_icsdT = df_MALM_icsd.drop('medianICSD',axis=1)
_ = surveyPRD.plot_PRD_effect_icsd(kmesh, 
                                    imin,
                                    df_MALM_icsdT, irr_log,
                                    ax=axs['b'],
                                    hours_interval=10,
                                    soil=False,
                                    # color=['red','blue'],
                                    detrend = 0.05
                                    )
# plt.legend(['SWC Left','SWC Right','ICSD Left','ICSD right'])
# ax_ERT = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,irr_log=irr_log,ax=ax2)
# surveyPRD.petro_plots(scalesurvey,k_indiv_merged,df_SWC,irr_log,ax=ax2)
# ax.set_title('icsd VS weight data')
axs['b'].set_ylabel('CSD')

ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

axs['b'].grid('on', which='major', axis='x',color='0.95' )
axs['b'].grid('on', which='major', axis='y',color='0.95' )
    
plt.savefig('../figures/icsd_PRD_effect.png', dpi=400)

#%%
# fig, ax = plt.subplots(1)

# scalesurvey.plot('datetime', 'weight_diff', style='.',ax=ax,
#                  color='k')
# # ax.set_ylim([-1e-3,1e-3])
# # weight_diff

# # dydx = np.gradient(y, dx)

# # dweight_dt = scalesurvey['weight (kg)'].diff() / scalesurvey['datetime'].index.to_series().diff().dt.total_seconds()
# dweight_dt = scalesurvey['weight (kg)'].diff() / scalesurvey['datetime'].diff().dt.total_seconds()


# dweight_dt.plot()

#%%
fig, ax = plt.subplots(1,sharex=True)
color = 'tab:red'
ax.set_xlabel('X-axis')
ax.set_ylabel('Y1-axis', color = color)
ax2 = ax.twinx()
color = 'tab:green'
ax2.set_ylabel('Y2-axis', color = color)     
# fig, axs = plt.subplots(3,sharex=False)
# ax_irr = surveyPRD.plot_timeline(ERT_log,irr_log,ax=axs[0],dropMALM=True)
ax_scale = surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax, color='red')
df_MALM_icsd = surveyPRD.plot_PRD_effect_icsd(kmesh, 
                                              imin,
                                              df_MALM_icsd, irr_log,
                                              ax=ax2,
                                              hours_interval=10)
    
# ax_ERT = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,irr_log=irr_log,ax=ax2)
# surveyPRD.petro_plots(scalesurvey,k_indiv_merged,df_SWC,irr_log,ax=ax2)
ax.set_title('icsd VS weight data')
plt.savefig('../figures/icsd_weight.png', dpi=400)


#%% DEPRECATED
# # normalized = (x-min(x))/(max(x)-min(x))
# # weightN_diff=((merge['weight (kg)']-merge['weight (kg)'].mean())/merge['weight (kg)'].std()).diff()
# # weightN_diff=merge['weight (kg)'].diff()
# # weightN_diff=(
# #                 (merge['weight (kg)']-merge['weight (kg)'].min())/
# #                 (merge['weight (kg)'].max()-merge['weight (kg)'].min())
# #             ).diff()

# # weightN_diff = merge['weight (kg)'].diff() / merge['datetime'].diff().dt.total_seconds()


# # meanResN_diff=((merge['meanRes']-merge['meanRes'].mean())/merge['meanRes'].std()).diff()
# meanResN_diff=((merge['meanRes']-merge['meanRes'].mean())/merge['meanRes'].std()).diff()
# # meanResN_diff=merge['meanRes'].diff()

# # weightN_diff = merge['meanRes'].diff() / merge['datetime'].diff().dt.total_seconds()



# # meanResN_diff=(
# #                 (merge['meanRes']-merge['meanRes'].min())/
# #                 (merge['weight (kg)'].max()-merge['weight (kg)'].min())
# #             ).diff()


# # meanResN_diff = (merge['meanRes']/max(merge['meanRes'])).diff()

# meanThetaN_diff=((merge['meanTheta']-merge['meanTheta'].mean())/merge['meanTheta'].std()).diff()

# # weightN_diff = (merge['weight (kg)']/max(merge['weight (kg)'])).diff()
# # meanResN_diff = (merge['meanRes']/max(merge['meanRes'])).diff()


# mosaic = 'ab'
# fig, axs = plt.subplot_mosaic(mosaic,
#                               constrained_layout=True,
#                               sharex=True
#                               )
# # fig, ax = plt.subplots(1,2)

# x11 = np.arange(-2,3)
# axs['a'].plot(x11,x11,
#          linestyle='--',
#          color='k')
# axs['a'].scatter(weightN_diff,
#             meanResN_diff,
#             color='k',
#             )
# axs['a'].set_xlabel(r'$\frac{d_{weight}}{dt}$', fontsize=14)
# axs['a'].set_ylabel(r'$\frac{d\rho}{dt}$', fontsize=14)
# axs['a'].axis('square')




# x11 = np.arange(-2,3)

# axs['b'].plot(x11,x11,
#          linestyle='--',
#          color='k')
# axs['b'].scatter(weightN_diff,
#             meanThetaN_diff,
#             color='k',
#             )
# axs['b'].set_xlabel(r'$\frac{d_{weight}}{dt}$', fontsize=14)
# axs['b'].set_ylabel(r'$\frac{d \Theta}{dt}$', fontsize=14)

# # axs['b'].set_ylabel('$d_{Theta}/dt$')
# axs['b'].axis('square')


# slope, intercept, r_value, p_value, std_err = stats.linregress(weightN_diff[1:], meanThetaN_diff[1:])
# yfit_theta = slope*weightN_diff[1:] + intercept
# R2_theta_text = '$R^{2}$=' + str(np.round(r_value,2))
# p_theta_text = '$p$=' + '{:.1e}'.format(p_value)


# # R2_res, p_res = proc.R2(weightN_diff[1:],meanResN_diff[1:])
# # R2_theta, p_theta = proc.R2(weightN_diff[1:],meanThetaN_diff[1:])

# from scipy.optimize import curve_fit
# def f(x, A, B): # this is your 'straight line' y=f(x)
#     return A*x + B

# # popt, pcov = curve_fit(f, weightN_diff[1:], meanResN_diff[1:]) # your data x, y to fit
# # yfit_res = popt[0]*weightN_diff[1:] + popt[1]


# slope, intercept, r_value, p_value, std_err = stats.linregress(weightN_diff[1:], meanResN_diff[1:])
# yfit_res = slope*weightN_diff[1:] + intercept
# R2_res_text = '$R^{2}$=' + str(np.round(r_value,2))
# p_res_text = '$p$=' + '{:.1e}'.format(p_value)

# axs['a'].plot(weightN_diff[1:],yfit_res,
#          linestyle='--',
#          color='y')

# # R2_res_fit, p_res_fit = proc.R2(weightN_diff[1:],yfit_res[1:])

# # R2_res_text = '$R^{2}$=' + str(np.round(R2_res,2))
# # p_res_text = '$p$=' + str(np.round(p_res,2))

# # R2_fit_res_text = '$R^{2}$=' + str(np.round(R2_res_fit,2))
# # p_fit_res_text = '$p$=' + str(np.round(p_res_fit,2))

# # R2_theta_text = '$R^{2}$=' + str(np.round(R2_theta,2))
# # p_theta_text = '$p$=' + str(np.round(p_theta,2))

# axs['a'].text(1,-2,R2_res_text)
# axs['a'].text(1,-2.2,p_res_text)
# # axs['a'].text(-2,1,R2_fit_res_text)
# # axs['a'].text(-2,0.8,p_fit_res_text)
# axs['b'].text(1,-2,R2_theta_text)
# axs['b'].text(1,-2.2,p_theta_text)


# for label, ax in axs.items():
#     # label physical distance in and down:
#     trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
#     ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
#             fontsize='medium', verticalalignment='top', fontfamily='serif',
#             bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))


# plt.savefig('../figures/weight_f_theta_Res_deprecated.png', dpi=400)


# # plt.scatter(merge['weight (kg)'],merge['meanRes'])

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
    
    

    