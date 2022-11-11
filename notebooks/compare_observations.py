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

import datetime
import surveyPRD

#%%

def get_cmd():
    parse = argparse.ArgumentParser()

    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[None], required=False) #[0,1,2,3,4,5]
    period.add_argument('-cycle', '--cycle', nargs='+',
                        help='list of cycle', type=int, default=[5,6,7,8], required=False)
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

#%% Load scale data + plot
# ----------------------------------------------------------------------------
scalesurvey = surveyPRD.load_scale_data(args.startD,args.endD)

fig, ax = plt.subplots(figsize=(10, 3),constrained_layout=True)
surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax)
ax.set_title('Scale raw data')
plt.savefig('../load_cells/figures/scale_raw_data.png', dpi=400)


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

# df_MALM_icsd = pd.read_csv(imaging_path + 
#                             inversionPathMALM + 
#                             'df_MALM_icsd' + 
#                             str([*args.cycle]) 
#                             + '.csv',
#                             header=[0, 1, 2],
#                             parse_dates=True,
#                             index_col=0
#                             )

# df_MALM_icsd = df_MALM_icsd.T.reset_index()
# df_MALM_icsd['datetime'] = df_MALM_icsd['datetime'].astype('datetime64[ns]')
# df_MALM_icsd.set_index(['filename', 'datetime', 'injection'], inplace=True)
# df_MALM_icsd = df_MALM_icsd.T


# imin = np.loadtxt(imaging_path + 
#                    inversionPathMALM + 
#                    'vrte_nodes_in_mesh.txt'
#                    )


#%% plot PRD_effect: ER
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
ax_ERT = surveyPRD.plot_PRD_effect_ER(kmesh,df_ERT,irr_log=irr_log,
                              unit='mS/m',ax=axs['b'])
for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
    
plt.savefig('../figures/Cond_variations.png', dpi=400)


#%% plot petro Archie
# -----------------------------------------------------------------------------

surveyPRD.petro_plots(
                        scalesurvey,
                        kmesh,
                        df_SWC,
                        irr_log
                      )
      
#%% plot ICSD
# -----------------------------------------------------------------------------

fig, ax = plt.subplots(1,sharex=True)
color = 'tab:red'
ax.set_xlabel('X-axis')
ax.set_ylabel('Y1-axis', color = color)
ax2 = ax.twinx()
color = 'tab:green'
ax2.set_ylabel('Y2-axis', color = color)     
# fig, axs = plt.subplots(3,sharex=False)
# ax_irr = surveyPRD.plot_timeline(ERT_log,irr_log,ax=axs[0],dropMALM=True)
# ax_scale = surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax, color='red')
ax_SWC = surveyPRD.plot_PRD_effect_SWC(kmesh,df_SWC,irr_log=irr_log,ax=ax2,
                                       LR=True, topdown=False,
                                       color=['green','yellow'])
df_MALM_icsd = surveyPRD.plot_PRD_effect_icsd(kmesh, 
                                              imin,
                                              df_MALM_icsd, irr_log,
                                              ax=ax2,
                                              hours_interval=10,
                                              soil=False,
                                              color=['red','blue'])
plt.legend(['SWC Left','SWC Right','ICSD Left','ICSD right'])
# ax_ERT = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,irr_log=irr_log,ax=ax2)
# surveyPRD.petro_plots(scalesurvey,k_indiv_merged,df_SWC,irr_log,ax=ax2)
ax.set_title('icsd VS weight data')
plt.savefig('../figures/icsd_swc.png', dpi=400)


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
    
    