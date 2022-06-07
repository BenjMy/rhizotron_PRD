#!/usr/bin/env python
# coding: utf-8

# In[1]:
from pyPRD import processing as proc
import numpy as np
import pyvista as pv
import argparse
import os 
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib as mpl
import math 

imaging_path = '../imaging_ERT_MALM/'


fmt = '%d/%m/%Y %H:%M'

# CYCLE 0
# -----------------------------
# startD = pd.to_datetime('13/05/2022 00:00', format=fmt)
# endD = pd.to_datetime('19/05/2022 17:00', format=fmt)


# CYCLE 1
# -----------------------------
# startD = pd.to_datetime('19/05/2022 16:47', format=fmt)
# endD = pd.to_datetime('25/05/2022 14:30', format=fmt)


# CYCLE 2
# -----------------------------
# startD = pd.to_datetime('25/05/2022 14:30', format=fmt)
# endD = pd.to_datetime('01/06/2022 15:50', format=fmt)


# CYCLE 4
# -----------------------------
# startD = pd.to_datetime('01/06/2022 01:00', format=fmt)
# endD = pd.to_datetime('02/06/2023 00:00', format=fmt)

# ALL CYCLES
# -----------------------------
startD = pd.to_datetime('13/05/2022 00:00', format=fmt)
endD = pd.to_datetime('02/06/2023 00:00', format=fmt)



#%%

reprocessed = True
def get_cmd():
    parse = argparse.ArgumentParser()
    
    period = parse.add_argument_group('period')
    period.add_argument('-period', type=str, help='period to analyse', default='2021', required=False)
    period.add_argument('-startD', type=str, help='start date to analyse', default=None, required=False)
    period.add_argument('-endD', type=str, help='end date', default=None, required=False)
    
    process_param = parse.add_argument_group('process_param')
    process_param.add_argument('-TL', type=bool, help='prepare Data Assimilation', default='True', required=False)
    process_param.add_argument('-irr', type=bool, help='run Data Assimilation', default='True', required=False)
    process_param.add_argument('-reprocessed', type=bool, help='reprocessed', default='True', required=False)
    reprocessed
    args = parse.parse_args()
    return(args)

#%%
import surveyPRD
ERT_log_all = surveyPRD.load_ERT_survey_log()


ERT_log = surveyPRD.load_ERT_survey_log(startDate=startD,endDate=endD)
irr_log = surveyPRD.load_irr_log_drive(startDate=startD,endDate=endD)


selected_files_ERT = ERT_log[ERT_log['method']=='ERT']['Name'].to_list()
selected_files_MALM = ERT_log[ERT_log['method']=='MALM']['Name'].to_list()

#%%
ax = surveyPRD.plot_timeline(ERT_log,irr_log)
ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=1))
ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))

    
#%%
# k.createBatchSurvey(imaging_path + 'filenames_ERT')
# k.invert(parallel=True)
    
k_indiv_merged = surveyPRD.process_ERT(imaging_path,
                                       selected_files_ERT,
                                       ERT_log,
                                       reprocessed=False,
                                       recip=5,
                                       )

#%%
# surveyPRD.add2pkl_ERT()
df_ERT = surveyPRD.add2df_ERT(k_indiv_merged,
                              selected_files_ERT,
                              ERT_log)

#%%
fig, ax = plt.subplots((2),figsize=(8,5), sharex=(True))
ax0 = surveyPRD.plot_timeline(ERT_log,irr_log, ax=ax[0], dropMALM=True)
# ax[0].set_ylim([-250,250])

# ax2 = ax.twinx()
# ax2.set_ylabel('Y2-axis')
# fig, ax = plt.subplots((2),figsize=(9,3), sharex=(True))
ax = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,ax=ax[1])

# # format xaxis with 4 month intervals
# ax.get_xaxis().set_major_locator(mdates.DayLocator(interval=1))
# #ax.get_xaxis().set_major_locator(mdates.MonthLocator(interval=0.1))
# ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m"))
# #plt.setp(ax.get_xticklabels(), rotation=30, ha="right")

    
plt.xticks(rotation=35)
plt.tight_layout()
plt.savefig(imaging_path+'inversionERT/'+ 'PRDeffect' + str([*selected_files_ERT]))

# # format xaxis with 4 month intervals
# ax.get_xaxis().set_major_locator(mdates.HourLocator(hours_interval))
# #ax.get_xaxis().set_major_locator(mdates.MonthLocator(interval=0.1))
# format_date='%d/%m/%Y %H:%M'
# ax.get_xaxis().set_major_formatter(mdates.DateFormatter(format_date))

#%%

# idTL = []
# for sf in selected_files_ERT:
#     idTL.append(ERT_log.index[ERT_log['Name']==sf].tolist()[0])

# # selected_files_ERT_reverse = list(reversed(selected_files_ERT))
# # idTL_reverse = list(reversed(idTL))

# # regType=1 # background constrainst inv

# k_TL = proc.invert_ERT_TL(
#                             imaging_path,
#                             files=selected_files_ERT,
#                             regType=1,
#                             recip=40,
#                             folderName=str([*idTL])
#                         )

# #%%
# proc.plot_ERT(k_TL[0], vmin=-10, vmax=1, attr="Sensitivity_map(log10)", index=1)
# proc.plot_ERT(k_TL[0], vmin=-10, vmax=1, attr="Sensitivity_map(log10)", index=0)

# proc.plot_ERT(k_TL[0], vmin=0, vmax=2, attr="Resistivity(log10)", index=0)
# proc.plot_ERT(k_TL[0], vmin=0, vmax=2, attr="Resistivity(log10)", index=1)
# proc.plot_ERT(k_TL[0], vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=0)
# proc.plot_ERT(k_TL[0], vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=1)
# proc.plot_ERT(k_TL[0], vmin=-20, vmax=20, attr="difference(percent)", index=1)

#%%
# df_ERT_TL = surveyPRD.add2df_ERT(k_TL,selected_files_ERT,ERT_log, TLflag=True)
# surveyPRD.plot_PRD_effect_ER(k_TL,df_ERT_TL)


# %%

# idt_irr = [[2, 3]]


# for itdIrr in idt_irr:
#     # selected_files_irr = [filenames_ERT[index] for index in itdIrr]
#     selected_files_irr = ERT_log_all[ERT_log_all['method']=='ERT'].iloc[itdIrr]['Name'].to_list()

#     k_TL = proc.invert_ERT_TL(
#         imaging_path,
#         files=selected_files_irr,
#         regType=1,
#         recip=10,
#         folderName= 'inversion_ERT_TL/BAirr_190522/'
#     )

#     proc.plot_ERT(k_TL, vmin=0, vmax=3, attr="Resistivity(log10)", index=0)
#     proc.plot_ERT(k_TL, vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=0)
#     proc.plot_ERT(k_TL, vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=1)
#     proc.plot_ERT(k_TL, vmin=-20, vmax=20, attr="difference(percent)", index=1)


#     k_TL = proc.invert_ERT_TL(
#         imaging_path,
#         files=selected_files_irr,
#         regType=2,
#         recip=2,
#         folderName= 'inversion_ERT_TL/BAirr_190522/'
#     )

#     proc.plot_ERT(k_TL, vmin=0, vmax=3, attr="Resistivity(log10)", index=0)
#     proc.plot_ERT(k_TL, vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=0)
#     proc.plot_ERT(k_TL, vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=1)
#     proc.plot_ERT(k_TL, vmin=-10, vmax=10, attr="difference(percent)", index=1)
#     k_TL.dirname


# df_ERT_TL_irr = surveyPRD.add2df_ERT(k_TL,selected_files_ERT,ERT_log)
# surveyPRD.plot_PRD_effect_ER(k_TL,df_ERT)



#%%
k_MALM = []
R_obs_stck = []
fig = plt.figure()
ax = fig.add_subplot()
for i, f in enumerate(selected_files_MALM):
    j = int(math.floor(i/2))    
    R_obs = proc.prepare_MALM(imaging_path,
                                f,
                                k_indiv_merged[j],
                                filter_seq=True,
                                percent_rec=1e9,
                                ax=ax
                                )
    R_obs_stck.append(R_obs)



#%%
proc.plot_scatter_MALM_diff(imaging_path,
                            R_obs_stck,
                            k_indiv_merged,
                            selected_files_MALM,
                            ERT_log
                            )
    
if len(selected_files_MALM) == 2*len(k_indiv_merged):
    proc.plot_scatter_MALM_diff_stem_soil(imaging_path,
                                     R_obs_stck,
                                     k_indiv_merged
                                     )

        
            
#%%
# k_indiv_merged.  


#     k_MALM.append(outMALM[0])
# nodes = outMALM[3]
# imin = outMALM[1]

# surveyPRD.add2df_MALM()

# %%

# m0 = []
# for i, f in enumerate(selected_files_MALM[0:2]):
#     j = ind_MALM_SOIL_STEM(i,selected_files_MALM,k_indiv_merged)
#     m0.append(proc.m0_MALM(imaging_path + 'inversionMALM/' + f,
#                            method_m0='F1')
#               )

#     # pl = proc.plot_ERT(k_indiv[0], vmin=0, vmax=3,
#     #                    attr="Resistivity(log10)",
#     #                    index=0, show=False)

#     proc.plot_m0_MALM(k_MALM[i], nodes, imin, m0[i])

# # %%

# sol = []
# for i, f in enumerate(selected_files_MALM[0:2]):
#     sol.append(proc.invert_MALM(imaging_path + 'inversionMALM/' + f,
#                                 wr=1)
#                )
#     proc.plot_MALM(k_MALM[i], nodes, imin, sol[i])

# # %%

# sol = []
# for i, f in enumerate(selected_files_MALM):
#     sol.append(proc.invert_MALM(imaging_path + 'inversionMALM/' + f,
#                                 pareto=True)
#                 )
#     proc.plot_MALM(k_MALM[i],nodes,imin,sol[i][0])



if __name__ == '__main__':

    '''
    - Read files
    - Invert ERT
    - Invert MALM
    '''
    
    # Read and parse args
    # -----------------------------------
    args = get_cmd()
    
    