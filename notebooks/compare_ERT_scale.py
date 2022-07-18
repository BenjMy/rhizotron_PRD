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


# surveyPRD.download_csv_scale()

#%%

def get_cmd():
    parse = argparse.ArgumentParser()
    
    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[4,5,6,7,8,9], required=False)
    period.add_argument('-cycle', type=list, help='cycle to analyse', default=[7], required=False)
    period.add_argument('-startD', type=str, help='start date to analyse', default=None, required=False)
    period.add_argument('-endD', type=str, help='end date', default=None, required=False)
    
    process_param = parse.add_argument_group('process_param')
    process_param.add_argument('-reprocessed', type=int, help='reprocessed', default=0, required=False)
    args = parse.parse_args()
    return(args)

args = get_cmd()
# args.cycle = [1,2]
args.cycle = [None]

#%%
imaging_path = '../imaging_ERT_MALM/'

# startD = pd.to_datetime('18/05/2022 00:00', format='%d/%m/%Y %H:%M')
# endD = pd.to_datetime('20/05/2022 00:00', format='%d/%m/%Y %H:%M')

fmt = '%d/%m/%Y,%H:%M'
# fmt = '%Y-%m-%d %H:%M:%S'
startD = None
endD = None


# startD = pd.to_datetime('01/06/2022 01:00', format='%d/%m/%Y %H:%M')
# endD = pd.to_datetime('02/06/2023 00:00', format='%d/%m/%Y %H:%M')


# # ALL CYCLES
# # -----------------------------
# startD = pd.to_datetime('13/05/2022 00:00', format=fmt)
# endD = pd.to_datetime('02/06/2023 00:00', format=fmt)



imaging_path = '../imaging_ERT_MALM/'

# %%

def get_cmd():
    parse = argparse.ArgumentParser()

    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[None], required=False) #[0,1,2,3,4,5]
    period.add_argument('-cycle', '--cycle', nargs='+',
                        help='list of cycle', type=int, default=[4,5,6,7,8,9], required=False)
    period.add_argument(
        '-startD', type=str, help='start date to analyse', default=None, required=False)
    period.add_argument('-endD', type=str, help='end date',
                        default=None, required=False)

    process_param = parse.add_argument_group('process_param')
    process_param.add_argument(
        '-recErr', type=int, help='Rec. error', default=5, required=False)
    process_param.add_argument(
        '-reprocessed', type=int, help='reprocessed', default=0, required=False)
    process_param.add_argument(
        '-TL', type=int, help='TimeLapse', default=0, required=False)
    process_param.add_argument(
        '-icsd', type=int, help='icsd proc.', default=1, required=False)
    process_param.add_argument(
        '-dim', type=str, help='icsd dim.', default='2d', required=False)
    process_param.add_argument(
        '-filter_seq', type=int, help='Filter sequence', default=1, required=False)
    process_param.add_argument('-filter_seq_rec', type=int,
                               help='Filter sequence rec', default=0, required=False)
    process_param.add_argument('-petro', type=int,
                               help='Apply Archie trans.', default=1, required=False)
    
    args = parse.parse_args()
    return(args)
#%%
args = get_cmd()
inversionPathMALM = surveyPRD.definePaths(args)
# ERTsurvey = surveyPRD.load_ERT_survey_log()
# d0 = ERT_log[ERT_log['Name']==selected_files_ERT[0]]['datetime'].dt.strftime(fmt)
# dend = ERT_log[ERT_log['Name']==selected_files_ERT[-1]]['datetime'].dt.strftime(fmt) 
args.startD = '08/06/2022,12:30'
args.endD = '12/07/2023,12:00'
# args.startD = '21/06/2022,12:00'
# args.endD = '30/06/2022,12:00'

ERT_log = surveyPRD.load_ERT_survey_log(startDate=args.startD,endDate=args.endD, 
                                        cycles=args.cycle)
irr_log = surveyPRD.load_irr_log_drive(startDate=args.startD,endDate=args.endD,
                                        cycles=args.cycle)
selected_files_ERT = ERT_log[ERT_log['method']=='ERT']['Name'].to_list()
selected_files_MALM = ERT_log[ERT_log['method']=='MALM']['Name'].to_list()

irr_log = surveyPRD.load_irr_log_drive(startDate=args.startD,endDate=args.endD)
scalesurvey = surveyPRD.load_scale_data(args.startD,args.endD)

#%%
# k.createBatchSurvey(imaging_path + 'filenames_ERT')
# k.invert(parallel=True)
    
k_indiv_merged = surveyPRD.process_ERT(imaging_path,
                                        selected_files_ERT,
                                        ERT_log,
                                        reprocessed=bool(args.reprocessed),
                                        recip=5,
                                        )
df_ERT = surveyPRD.add2df_ERT(k_indiv_merged,selected_files_ERT,ERT_log)
df_ERT_res = (1/df_ERT)*1e3
df_ERT_res.columns
# df_ERT


#%%
fig, ax = plt.subplots(figsize=(10, 3),constrained_layout=True)
surveyPRD.plot_weight_dynamic(scalesurvey,ax=ax)
ax.set_title('Scale raw data')
plt.savefig('../load_cells/figures/scale_raw_data.png', dpi=400)

#%%
ax = surveyPRD.plot_timeline(ERT_log,irr_log,
                        dropMALM=True)
ax.set_title('Timeline irr. pattern')
plt.savefig('../figures/irr_pattern.png', dpi=400)

#%%

mosaic = [
            ['a.'],['b.']
          ]
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=True)

# fig, ax = plt.subplots(2,1,figsize=(7, 3),constrained_layout=True,
#                        sharex=True)
# axs['a.'].set_title()
surveyPRD.plot_timeline(ERT_log,irr_log,
                        dropMALM=True,ax=axs['a.'])
surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,irr_log=irr_log,
                             unit='mS/m',ax=axs['b.'])
# ax[1].set_title('Cond. variations')

for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
    
plt.savefig('../figures/Cond_variations.png', dpi=400)


#%%
fig, axs = plt.subplot_mosaic(mosaic,
                              constrained_layout=True,
                              sharex=True)
surveyPRD.plot_timeline(ERT_log,irr_log,
                        dropMALM=True,ax=axs['a.'])
surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT_res,irr_log, unit='Ohm.m',
                             ax=axs['b.'])

for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))
    
ax.set_title('ER variations')
plt.savefig('../figures/ER_variations.png', dpi=400)


#%%
# fig, axs = plt.subplot_mosaic(mosaic,
#                               constrained_layout=True,
#                               sharex=True)
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
ax_ERT = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,irr_log=irr_log,ax=ax2)
ax.set_title('ER VS weight data')
plt.savefig('../figures/ER_weight.png', dpi=400)

#%%

if args.petro:
    rFluid = 1
    porosity = 0.4
    df_SWC = df_ERT_res.copy()
    for cols in list(df_ERT.columns):
        if type(cols)==str:
            pass
        else:
            df_SWC[cols] = surveyPRD.Archie_rho2sat(df_ERT_res[cols].to_numpy(),
                                                      rFluid, porosity, 
                                                      a=1.0, m=2.0, n=2.0)
            
surveyPRD.petro_plots(scalesurvey,k_indiv_merged,df_SWC,irr_log)
      
# %% ICSD analysis

args.cycle=[4,5,6,7,8,9]
df_MALM_icsd = pd.read_csv(imaging_path + 
                            inversionPathMALM + 
                            'df_MALM_icsd' + 
                            str([*args.cycle]) 
                            + '.csv',
                            header=[0, 1])

imin = np.loadtxt(imaging_path + 
                   inversionPathMALM + 
                   'vrte_nodes_in_mesh.txt'
                   )


# for tc in zip(df_MALM_icsd.columns):
#     print(tc)
#     if 'LR' in tc[0]:
#         continue
#         if tc[0][2] is not None:
#             print(tc[0][2])
        
            
# df_MALM_icsd = pd.read_csv('df_MALM_icsd3456789'
#                             + '.csv')

# df_SWC

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
df_MALM_icsd = surveyPRD.plot_PRD_effect_icsd(k_indiv_merged, 
                                              imin,
                                              df_MALM_icsd, irr_log,
                                              ax=ax,
                                              hours_interval=10)
    
ax_ERT = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,irr_log=irr_log,ax=ax2)
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
    
    