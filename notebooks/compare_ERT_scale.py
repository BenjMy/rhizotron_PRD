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

import surveyPRD
#%%
imaging_path = '../imaging_ERT_MALM/'

# startD = pd.to_datetime('18/05/2022 00:00', format='%d/%m/%Y %H:%M')
# endD = pd.to_datetime('20/05/2022 00:00', format='%d/%m/%Y %H:%M')

fmt = '%d/%m/%Y %H:%M'

startD = None
endD = None


startD = pd.to_datetime('01/06/2022 01:00', format='%d/%m/%Y %H:%M')
endD = pd.to_datetime('02/06/2023 00:00', format='%d/%m/%Y %H:%M')


# ALL CYCLES
# -----------------------------
startD = pd.to_datetime('13/05/2022 00:00', format=fmt)
endD = pd.to_datetime('02/06/2023 00:00', format=fmt)


ERTsurvey = surveyPRD.load_ERT_survey_log()
ERT_log = surveyPRD.load_ERT_survey_log(startDate=startD,endDate=endD)
selected_files_ERT = ERT_log[ERT_log['method']=='ERT']['Name'].to_list()
selected_files_MALM = ERT_log[ERT_log['method']=='MALM']['Name'].to_list()


irr_log = surveyPRD.load_irr_log_drive(startDate=startD,endDate=endD)
scalesurvey = surveyPRD.load_scale_data(startD,endD)


#%%
# k.createBatchSurvey(imaging_path + 'filenames_ERT')
# k.invert(parallel=True)
    
k_indiv_merged = surveyPRD.process_ERT(imaging_path,
                                        selected_files_ERT,
                                        ERT_log,
                                        reprocessed=False,
                                        recip=5,
                                        )
# # surveyPRD.add2pkl_ERT()
df_ERT = surveyPRD.add2df_ERT(k_indiv_merged,selected_files_ERT,ERT_log)

#%%
surveyPRD.plot_weight_dynamic(scalesurvey)


#%%
surveyPRD.plot_timeline(ERT_log,irr_log)

#%%
surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT)

#%%
fig, axs = plt.subplots(3,sharex=True)
# fig, axs = plt.subplots(3,sharex=False)

ax_irr = surveyPRD.plot_timeline(ERT_log,irr_log,ax=axs[0])
ax_scale = surveyPRD.plot_weight_dynamic(scalesurvey,ax=axs[1])
ax_ERT = surveyPRD.plot_PRD_effect_ER(k_indiv_merged,df_ERT,ax=axs[2])
