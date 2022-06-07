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
import matplotlib.dates as mdates
from resipy import Project


def plot_PRD_effect_ER(k_indiv_merged,df_ERT,ax=None,
                       hours_interval=10):
    
    if ax == None:
        fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)
        
    meshXYZ = k_indiv_merged[0].mesh.df[['X','Y','Z']]
    meshXYZ['X'].min()
    meshXYZ['X'].max()
    
    idx_mesh_left = meshXYZ['X']<meshXYZ['X'].max()/2
    idx_mesh_right = meshXYZ['X']>meshXYZ['X'].max()/2
    df_ERT['LR'] = None
    df_ERT['LR'][idx_mesh_left]= 'Left'
    df_ERT['LR'][idx_mesh_right]= 'Right'
    
    df_groups = df_ERT.groupby(['LR']).mean()
    # df_groups.T.plot(xlabel='date',ylabel='mean Conductivity (mS/m)', ax=ax)
    # df_groups.T.plot.scatter(xlabel='date',ylabel='mean Conductivity (mS/m)', ax=ax)
    # df_groups.T.plot(xlabel='date',ylabel='mean Conductivity (mS/m)', ax=ax, kind='bar')
    df_groups.T.plot(xlabel='date',ylabel='mean Conductivity (mS/m)', ax=ax)


    return ax


def plot_timeline(ERT_log,irr_log, ax=None,
                  show=False, **kwargs):
    dates_ERT = ERT_log['datetime']
    names = ERT_log['method']
    
    if 'dropMALM' in kwargs:
        dates_ERT.drop(ERT_log[ERT_log['method']=='MALM'].index.to_list(), axis=0, inplace=True)
        names = ['ERT + MALM']*len(dates_ERT)
    
    dates_irr = irr_log['datetime']
    ml_irr = irr_log['quantity (mL)']
        
        
    id_left = irr_log[irr_log['where'].str.contains("1")].index.tolist()
    id_right = irr_log[irr_log['where'].str.contains("8")].index.tolist()
    
    ml_irr_signed = ml_irr
    ml_irr_signed[id_left] = ml_irr[id_left]*-1
    
    # Choose some nice levels
    # levels = np.tile([-5, 5, -3, 3, -1, 1],
    #                  int(np.ceil(len(dates)/6)))[:len(dates)]
    
    
    # Create figure and plot a stem plot with the date
    if ax == None:
        fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)
    ax.set(title="Timeline rhizotron experiment")
    
    levels_irr = np.tile(ml_irr_signed.to_list(),
                     int(np.ceil(len(dates_irr))))[:len(dates_irr)]
    levels_ERT = np.tile([1],
                     int(np.ceil(len(dates_ERT))))[:len(dates_ERT)]
    
    markerline_irr, stemline_irr, baseline_irr = ax.stem(dates_irr, levels_irr,
                                             linefmt="--", basefmt="k-")
    markerline_ERT, stemline_ERT, baseline_ERT = ax.stem(dates_ERT, levels_ERT,
                                             linefmt="C4-", basefmt="k-")
    
    # plt.setp(markerline, mec="k", mfc="w", zorder=3)
    
    # Shift the markers to the baseline by replacing the y-data by zeros.
    #markerline.set_ydata(np.zeros(len(dates)))
    
    # annotate lines
    vert = np.array(['top', 'bottom'])[(levels_ERT > 0).astype(int)]
    for d, l, r, va in zip(dates_ERT, levels_ERT, names, vert):
        ax.annotate(r, xy=(d, l), xytext=(-3, np.sign(l)*3),
                    textcoords="offset points", va=va, ha="right",rotation=65)
       
    # # remove y axis and spines
    # ax.get_yaxis().set_visible(False)
    # for spine in ["left", "top", "right"]:
    #     ax.spines[spine].set_visible(False)
    
    # ax.margins(y=0.1)
    if show:
        plt.show()
    return ax
    
def select_from_dates(df,start,end):
    mask = (df['datetime'] > start) & (df['datetime'] <= end)
    filtered_df=df.loc[mask]
    return filtered_df

#%%
# Load ERT data filenames

def add2pkl_ERT():
    pass

def add2df_ERT(kresipy,selected_files_ERT,ERT_log, TL_flag=False):
    ''' Create a dataframe containing datetimes and the respectives ER values after inversion'''
    df_ERT = []
    for i, f in enumerate(selected_files_ERT):
        date_f = ERT_log['datetime'].loc[ERT_log[ERT_log['Name']==f].index.to_list()[0]]
        kresipy[i].meshResults[0].df['Conductivity(mS/m)']
        if len(df_ERT)==0:
            df_ERT = pd.DataFrame({date_f:kresipy[i].meshResults[0].df['Conductivity(mS/m)']})
        else:
            df_ERT2merge = pd.DataFrame({date_f:kresipy[i].meshResults[0].df['Conductivity(mS/m)']})
            df_ERT = pd.concat([df_ERT,df_ERT2merge],axis=1)
    return df_ERT

def load_ERT_survey_log_drive():
    sheet_id = '1poMI4lgeolZDPdScZyoUIDDZfGHUrsyBgbiyCU0HRI8'
    # sheet_name = ['2nd_run']
    # data = []
    # for sn in sheet_name:
    #     url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sn}'
    #     data = pd.read_csv(url, decimal=',',header='infer', nrows=1)
    #     print(url)
    pass


def load_ERT_survey_log(csv2read=('/home/ben/Documents/GitHub/BenjMy/'+
                                  'rhizotron_PRD/imaging_ERT_MALM/'+
                                  'PRD_Measurements_log - 2nd_run.csv'),
                        startDate=None,
                        endDate=None,
                        ):
    survey_log = pd.read_csv(csv2read, 
                             decimal=',',
                             header='infer'
                             )
    survey_log=survey_log.dropna(subset=['Instrument','Date'])
    survey_log['Date']
    survey_log['Time (CET in hh:mm)']
    
    survey_log['datetime_tmp'] = survey_log['Date'] + ' ' + survey_log['Time (CET in hh:mm)']
    # survey_log['datetime'] = pd.to_datetime(survey_log['datetime'], format='%d/%-m/%Y %H:%M')
    survey_log['datetime'] = pd.to_datetime(survey_log['datetime_tmp'], format='%d/%m/%Y %H:%M')
    
    id_ERT = survey_log[survey_log['Name'].str.contains("ERT")].index.tolist()
    id_MALM = survey_log[survey_log['Name'].str.contains("MALM")].index.tolist()

    # survey_log['method'] = 
    survey_log["method"] = None
    survey_log['method'][id_ERT] = 'ERT'
    survey_log['method'][id_MALM] = 'MALM'

    if startDate is not None:
        survey_log = select_from_dates(survey_log,startDate,endDate)

    return survey_log


def add_process_log(df,f):
    df['inverted_bool'] = np.ones(len(df))*False
    df[df['Name']==f]['inverted_bool'] = True
    return df

def indiv_invert_ERT(imaging_path,f,recip=5):
    k_indiv = proc.invert_ERT(imaging_path,
                                    filename=f,
                                    recip=recip
                                )
    proc.plot_ERT(k_indiv, vmin=0, vmax=50,
                  attr="Resistivity(ohm.m)",
                  index=0)
    return k_indiv

def process_ERT(imaging_path,selected_files_ERT,ERT_log,reprocessed=False,
                recip=5):
    k_indiv_merged=[]
    for i, f in enumerate(selected_files_ERT):
        print(f)
        date_f = ERT_log['datetime'].loc[ERT_log[ERT_log['Name']==f].index.to_list()[0]]
        if os.path.exists(os.path.join(imaging_path,'inversionERT') + '/'+ f + "bkp.resipy"):
            if reprocessed:
               k_indiv = indiv_invert_ERT(imaging_path,f,recip=recip)
               k_indiv_merged.append(k_indiv)
            else:
               print('reload file from *.resipy backup')
               k_loaded = Project(dirname=imaging_path, typ='R3t')
               k_loaded.loadProject(os.path.join(imaging_path,'inversionERT') + '/'+ f + "bkp.resipy")
               k_indiv_merged.append(k_loaded)
               pass
        else:
            k_indiv = indiv_invert_ERT(imaging_path,f,recip=recip)
            k_indiv_merged.append(k_indiv)
    return k_indiv_merged


def add2df_MALM():
    pass


#%%


def load_irr_log_drive(startDate=None,
                       endDate=None,
                       ):
    # sheet_id = '15AvYhQywK04VDfoBbAnqXP_YJeBmIwi-AxYzuRNjYdo'
    # sheet_name = ['2nd_run']
    # data = []
    # for sn in sheet_name:
    #     url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sn}'
    #     data = pd.read_csv(url)
    #     print(data)
        
    # sheet_url = "https://docs.google.com/spreadsheets/d/{sheet_id}/edit#gid=0"
    # url_1 = sheet_url.replace("/edit#gid=", "/export?format=csv&gid=")
    csv2read = '/home/ben/Documents/GitHub/BenjMy/rhizotron_PRD/Irrigation_log - 2nd_run.csv'
    irr_log = pd.read_csv(csv2read, decimal=',')  

    irr_log['datetime_tmp'] = irr_log['date'] + ' ' + irr_log['start']
    irr_log['datetime'] = pd.to_datetime(irr_log['datetime_tmp'], format='%d/%m/%Y %H:%M')

    if startDate is not None:
        irr_log = select_from_dates(irr_log,startDate,endDate)

    return irr_log

# Load scale data filenames

def plot_weight_dynamic(scalesurvey,ax=None):
    
    if ax == None:
        fig, ax = plt.subplots(figsize=(18, 5))
    
    scalesurvey.plot('datetime', 'weight (kg)', style='o', ax=ax)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
    #ax.xaxis.set_minor_locator(mdates.HourLocator(interval=4))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
    #ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))

    # Set title and labels for axes
    ax.set(xlabel="Date",
           ylabel="weight (kg)",
           title="Raw data sheet")
    return ax

def load_scale_data(startDate=None,
                    endDate=None,
                    ):

    sheet_id = '1GqDKU8fC0Pw_DQpUD-MBLpN-AXWrAzLuO1hABxxYIhc'
    sheet_name = ['apr_6-11','apr_11-13','apr_13-21','apr_21-26','apr_28-may_2',
                  'may_2-9','may_9-12','may_13-16','may_16-19','may_19-23',
                  'may_23-31',
                  ]

    list_data_scale = []
    for sn in sheet_name:
        url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sn}'
        data_scale_i = pd.read_csv(url, decimal=',')
        
        data_scale_i['abs_date'] = pd.to_datetime(data_scale_i['date'] + " " + data_scale_i['time'],
                                            infer_datetime_format=True)  
        datestr = data_scale_i['date'] + " " + data_scale_i['time']
        data_scale_i['abs_date'] = pd.to_datetime(data_scale_i['date'] + " " + data_scale_i['time'],
                                        format="%d/%m/%Y %H:%M")         
        initial_date = data_scale_i['abs_date'][0]       
        dates = []
        for i in range(len(data_scale_i['time'])):
            dates.append(initial_date + pd.Timedelta(seconds=data_scale_i['sec'][i]))
        date_time = [t.strftime("%Y-%m-%d %H:%M:%S") for t in dates]
        data_scale_i['datetime']=date_time
        data_scale_i['datetime']= pd.to_datetime(data_scale_i['datetime'])
        list_data_scale.append(data_scale_i)
        


    df = pd.concat(list_data_scale, axis=0, ignore_index=False)
    
    if startDate is not None:
        df = select_from_dates(df,startDate,endDate)


    return df



# Plot irrigation schedule (see Legnaro plot)



if __name__ == '__main__':

    '''
    - Read files
    - Invert ERT
    - Invert MALM
    '''
    
    # Read and parse args
    # -----------------------------------
    args = get_cmd()
    
    