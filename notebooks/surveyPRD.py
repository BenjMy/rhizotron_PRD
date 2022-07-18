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
import math



def definePaths(args):
    
    inversionPathMALM = 'inversionMALM/'
    if args.filter_seq:
        inversionPathMALM = 'inversionMALM_filter/'
    if args.filter_seq_rec:
        inversionPathMALM = 'inversionMALM_filter_rec/'
        
    return inversionPathMALM


def load_log_files(args):
    ERT_log = load_ERT_survey_log(startDate=args.startD, endDate=args.endD,
                                  cycles=args.cycle)
    irr_log = load_irr_log_drive(startDate=args.startD, endDate=args.endD,
                                 cycles=args.cycle)

    selected_files_ERT = ERT_log[ERT_log['method'] == 'ERT']['Name'].to_list()
    selected_files_MALM = ERT_log[ERT_log['method']
                                  == 'MALM']['Name'].to_list()

    # %%
    ax = plot_timeline(ERT_log, irr_log, dropMALM=True)
    tdelta = (ERT_log['datetime'].iloc[0]-ERT_log['datetime'].iloc[-1])/10
    ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=55))
    ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))

    return selected_files_ERT, selected_files_MALM, ERT_log, irr_log


def set_Xlim_LR_from_mesh(k_indiv_merged):
    ''' Define lim from left to right using mesh nodes positions'''
    meshXYZ = k_indiv_merged[0].mesh.df[['X', 'Y', 'Z']]
    meshXYZ['X'].min()
    meshXYZ['X'].max()
    return meshXYZ, meshXYZ['X'].max()/2


def plot_PRD_effect_icsd(k_indiv_merged, vrte_in_mesh, 
                         df_MALM_icsd, irr_log,
                         ax=None,
                         hours_interval=10,
                         **kwargs):

        
    ylabel = r'Normalised Current Density'
    if 'ylabel' in kwargs:
        ylabel = kwargs['ylabel']
        
        
    if ax == None:
        fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)

    id_left = irr_log[irr_log['where'].str.contains("1")].index.tolist()
    id_right = irr_log[irr_log['where'].str.contains("8")].index.tolist()

    meshXYZ, midlleX = set_Xlim_LR_from_mesh(k_indiv_merged)

    idx_vrte_left_bool = meshXYZ['X'][vrte_in_mesh] < midlleX
    idx_vrte_left_bool = idx_vrte_left_bool.reset_index()
    idx_vrte_left = idx_vrte_left_bool[idx_vrte_left_bool['X']].index.values
        
    # df_MALM_icsd = df_MALM_icsd.drop(['LR'], axis=1)
    # df_MALM_icsd[('LR',None,None)] = 'Right'
    
    
    left_right = ['Right']*len(df_MALM_icsd)
    for idxl in idx_vrte_left:
        left_right[idxl] = 'Left'
    
    df_MALM_icsd['LR']= left_right
    # df_MALM_icsd[('LR',None,None)] = 'Right'

    # df_groups = df_MALM_icsd.groupby(['LR'],level=0).mean()
    # df_groups.T.plot(xlabel='date', ylabel='mean Conductivity (mS/m)', ax=ax,
    #                  linestyle='--', marker='+', color=['green', 'orange'])
    
    tuple_cols_2_select_soil = []
    tuple_cols_2_select_stem = []
    
    for tc in zip(df_MALM_icsd.columns):
        print(tc[0][2] )
        # if 'LR' in tc[0]:
        #     continue
        if tc[0][2] is not None:
            if ('Soil' in tc[0][2]):
                tuple_cols_2_select_soil.append(tc[0])
            if ('Stem' in tc[0][2]):
                tuple_cols_2_select_stem.append(tc[0])

    
    if len(left_right) == len(df_MALM_icsd):

        soil_inject_df_MALM = df_MALM_icsd[tuple_cols_2_select_soil]
        soil_inject_df_MALM['LR']= left_right
        df_groups_soil_inj = soil_inject_df_MALM.droplevel(level=[1,2],axis=1).groupby(['LR']).sum()
        df_groups_soil_inj.columns = df_MALM_icsd[tuple_cols_2_select_soil].columns

        df_groups_soil_inj.T.droplevel(level=[0,2]).plot(xlabel='date', ylabel=ylabel, ax=ax,
                                linestyle=':',
                                marker='o',
                                label='soil inj. sum curr.',
                                color=['grey', 'grey']
                                )


        stem_inject_df_MALM = df_MALM_icsd[tuple_cols_2_select_stem]
        stem_inject_df_MALM['LR']= left_right

        df_groups_stem_inj = stem_inject_df_MALM.droplevel(level=[1,2],axis=1).groupby(['LR']).sum()
        # print(stem_inject_df_MALM)
        df_groups_stem_inj.columns = df_MALM_icsd[tuple_cols_2_select_stem].columns

        df_groups_stem_inj.T.droplevel(level=[0,2]).plot(xlabel='date', ylabel=ylabel, ax=ax,
                                    linestyle='--',
                                    marker='v',
                                    label='stem inj. sum curr.',
                                    color=['orange', 'green'])
            
    return ax

def plot_PRD_effect_SWC(k_indiv_merged, df_SWC, irr_log,
                       ax=None,
                       hours_interval=10,
                       unit='m3/m3',
                       **kwargs,
                       ):
    
    ax = plot_PRD_effect(k_indiv_merged, df_SWC, irr_log,
                           ax,
                           hours_interval,
                           unit=unit,
                           **kwargs,
                         )
    return ax
    
def plot_PRD_effect_ER(k_indiv_merged, df_ERT, irr_log,
                       ax=None,
                       hours_interval=10,
                       unit='mS/m'):
    
    ax = plot_PRD_effect(k_indiv_merged, df_ERT, irr_log,
                           ax,
                           hours_interval,
                           unit=unit)
    return ax

def plot_PRD_effect(k_indiv_merged, df, irr_log,
                       ax=None,
                       hours_interval=10,
                       unit='mS/m',
                       topdown=False,
                       **kwargs):

    ylabel = 'mean Conductivity (mS/m)'
    if unit == 'Ohm.m':
        ylabel = 'mean Resistivity (Ohm.m)'
    elif unit == 'm3/m3':
        ylabel = 'mean SWC (m3/m3)'

    LR = True
    if 'LR' in kwargs:
        LR = kwargs['LR']
    if 'topdown' in kwargs:
        topdown = kwargs['topdown']

            
    if ax == None:
        fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)

    id_left = irr_log[irr_log['where'].str.contains("1")].index.tolist()
    id_right = irr_log[irr_log['where'].str.contains("8")].index.tolist()

    meshXYZ, midlleX = set_Xlim_LR_from_mesh(k_indiv_merged)
    midlleZ = meshXYZ['Z'].max()/2
    
    idx_mesh_left = meshXYZ['X'] < midlleX
    idx_mesh_right = meshXYZ['X'] > midlleX
    df['LR'] = None
    df['LR'][idx_mesh_left] = 'Left'
    df['LR'][idx_mesh_right] = 'Right'


    idx_mesh_top = meshXYZ['Y'] < midlleZ
    idx_mesh_down = meshXYZ['Z'] > midlleZ
    df['TopDown'] = None
    df['TopDown'][idx_mesh_top] = 'Top'
    df['TopDown'][idx_mesh_down] = 'Down'
    
    
    # df_ERT['color_plot'] = None
    # df_ERT['color_plot'][idx_mesh_left]= 'green'
    # df_ERT['color_plot'][idx_mesh_right]= 'orange'

    # df_ERT['color_plot'] = 'orange'

    # irr_log['color_plot'] = 'orange'
    # try:
    #     for idl in id_left:
    #         irr_log['color_plot'][idl] = 'green'
    # except:
    #     pass

    if LR and topdown:
        df_groups = df.groupby(['LR','TopDown']).mean()
        df_groups.T.plot(xlabel='date', ylabel=ylabel, ax=ax,
                         linestyle='-.', marker='+', color=['limegreen', 'darkgreen','bisque', 'darkorange'])
    elif LR:
        df_groups = df.groupby(['LR']).mean()
        df_groups.T.plot(xlabel='date', ylabel=ylabel, ax=ax,
                         linestyle='-.', marker='+', color=['darkgreen','darkorange'])
    elif topdown:
        df_groups = df.groupby(['TopDown']).mean()
        df_groups.T.plot(xlabel='date', ylabel=ylabel, ax=ax,
                         linestyle='-.', marker='+', color=['darkgreen','darkorange'])
    else:
        df.T.plot(xlabel='date', ylabel=ylabel, ax=ax,
                         linestyle='-.', marker='+', color=['darkgreen'])
        
        
    # plt.grid(axis='x', color='0.95')
    ax.grid('on', which='major', axis='x',color='0.95' )
    ax.grid('on', which='major', axis='y',color='0.95' )


    # ax.fill_between(df_groups.columns,     
    #                 min(df_groups.min(axis=0)),
    #                 max(df_groups.max(axis=0)),
    #                 where=df_groups[df_groups.columns[0]].values>1,
    #                 color='green', alpha=0.5, transform=ax.get_xaxis_transform())

    # ax.fill_between(df_groups.columns, 0, 1, where=y > threshold,
    #                 color='white', alpha=0.5, transform=ax.get_xaxis_transform())

    return ax


def plot_PRD_effect_MALM_diff(k_indiv_merged, df_MALM, ax=None,
                              hours_interval=10, **kwargs):

    print(df_MALM)

    if ax == None:
        fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)

    ylabel = r'Diff R ($\Omega$)'
    if 'ylabel' in kwargs:
        ylabel = kwargs['ylabel']

    tuple_cols_2_select_soil_diff = []
    tuple_cols_2_select_stem_diff = []
    for tc in zip(df_MALM.columns):
        print(tc)
        if tc[0][2] is not False:
            if ('diff_soil' in tc[0][2]):
                tuple_cols_2_select_soil_diff.append(tc[0])
            if ('diff_stem' in tc[0][2]):
                tuple_cols_2_select_stem_diff.append(tc[0])

    soil_inject_df_MALM_diff = df_MALM[tuple_cols_2_select_soil_diff]
    soil_inject_df_MALM_diff['LR'] = df_MALM['LR']

    df_groups_soil_inj = soil_inject_df_MALM_diff.groupby(['LR']).mean()
    df_groups_soil_inj.T.droplevel(level=[1, 2]).plot(xlabel='date', ylabel=ylabel, ax=ax,
                                                      linestyle=':',
                                                      marker='o',
                                                      label='soil inj.',
                                                      color=['grey', 'grey'])

    stem_inject_df_MALM_diff = df_MALM[tuple_cols_2_select_stem_diff]
    stem_inject_df_MALM_diff['LR'] = df_MALM['LR']

    df_groups_stem_inj = stem_inject_df_MALM_diff.groupby(['LR']).mean()
    df_groups_stem_inj.T.droplevel(level=[1, 2]).plot(xlabel='date', ylabel=ylabel, ax=ax,
                                                      linestyle='--',
                                                      marker='v',
                                                      label='stem inj.',
                                                      color=['orange', 'green'])

    return ax

    # if ax == None:
    #     fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)
    # _, midlleX = set_Xlim_LR_from_mesh(k_indiv_merged)

    # ylabel=r'R ($\Omega$)'
    # if 'ylabel' in kwargs:
    #     ylabel = kwargs['ylabel']

    # elecsMALM = k_indiv_merged[0].surveys[0].elec
    # elecsMALM = elecsMALM.drop(index=[64,70,71], axis=0)
    # xelecs = elecsMALM['x']

    # idx_elecs_left = xelecs<midlleX
    # idx_elecs_right = xelecs>midlleX

    # LR = ['Left']*len(df_MALM)
    # for i, right in enumerate(idx_elecs_right):
    #     if right:
    #         LR[i] = 'Right'
    # df_MALM['LR'] = LR


def plot_PRD_effect_MALM(k_indiv_merged, df_MALM, ax=None,
                         hours_interval=10, **kwargs):

    if ax == None:
        fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)
    _, midlleX = set_Xlim_LR_from_mesh(k_indiv_merged)

    ylabel = r'R ($\Omega$)'
    if 'ylabel' in kwargs:
        ylabel = kwargs['ylabel']

    elecsMALM = k_indiv_merged[0].surveys[0].elec
    elecsMALM = elecsMALM.drop(index=[64, 70, 71], axis=0)
    xelecs = elecsMALM['x']

    idx_elecs_left = xelecs < midlleX
    idx_elecs_right = xelecs > midlleX

    LR = ['Left']*len(df_MALM)
    for i, right in enumerate(idx_elecs_right):
        if right:
            LR[i] = 'Right'
    df_MALM['LR'] = LR

    if len(idx_elecs_left) == len(df_MALM):
        # if df_MALM.columns.ndims
        try:
            tuple_cols_2_select_soil = []
            tuple_cols_2_select_stem = []
            for tc in zip(df_MALM.columns):
                print(tc[0][0])
                if ('Soil' in tc[0][1]) | (type(tc[0][0]) == str):
                    tuple_cols_2_select_soil.append(tc[0])
                if ('Stem' in tc[0][1]) | (type(tc[0][0]) == str):
                    tuple_cols_2_select_stem.append(tc[0])

            soil_inject_df_MALM = df_MALM[tuple_cols_2_select_soil]
            df_groups_soil_inj = soil_inject_df_MALM.groupby(['LR']).mean()
            df_groups_soil_inj.T.droplevel(level=1).plot(xlabel='date', ylabel=ylabel, ax=ax,
                                                         linestyle=':',
                                                         marker='o',
                                                         label='soil inj.',
                                                         color=['grey', 'grey'])

            stem_inject_df_MALM = df_MALM[tuple_cols_2_select_stem]
            df_groups_stem_inj = stem_inject_df_MALM.groupby(['LR']).mean()
            df_groups_stem_inj.T.droplevel(level=1).plot(xlabel='date', ylabel=ylabel, ax=ax,
                                                         linestyle='--',
                                                         marker='v',
                                                         label='stem inj.',
                                                         color=['orange', 'green'])
        except:
            df_MALM_groupLR = df_MALM.groupby(['LR']).mean()
            df_MALM_groupLR.T.plot(xlabel='date', ylabel=ylabel, ax=ax,
                                   linestyle='--',
                                   marker='v',
                                   label='stem inj.',
                                   color=['orange', 'green'])

    else:
        # raise ValueError('Impossible to split the seq')
        print('Impossible to split the seq')

    return ax


def plot_timeline(ERT_log, irr_log, ax=None,
                  show=False, **kwargs):
    dates_ERT = ERT_log['datetime']
    names = ERT_log['method']

    if 'dropMALM' in kwargs:
        dates_ERT.drop(ERT_log[ERT_log['method'] ==
                       'MALM'].index.to_list(), axis=0, inplace=True)
        # names = ['ERT + MALM']*len(dates_ERT)
        names = list(np.arange(0, len(dates_ERT)))
        names_str = list(map(str, names))

    dates_irr = irr_log['datetime']
    ml_irr = irr_log['quantity (mL)']

    id_left = irr_log[irr_log['where'].str.contains("1")].index.tolist()
    id_right = irr_log[irr_log['where'].str.contains("8")].index.tolist()

    ml_irr_signed = ml_irr
    # ml_irr_signed[id_left] = ml_irr[id_left]*-1

    # Choose some nice levels
    # levels = np.tile([-5, 5, -3, 3, -1, 1],
    #                  int(np.ceil(len(dates)/6)))[:len(dates)]

    # Create figure and plot a stem plot with the date
    if ax == None:
        fig, ax = plt.subplots(figsize=(8.8, 4), constrained_layout=True)
    # ax.set(title="Timeline rhizotron experiment")

    levels_irr = np.tile(ml_irr_signed.to_list(),
                         int(np.ceil(len(dates_irr))))[:len(dates_irr)]
    levels_ERT = np.tile([1],
                         int(np.ceil(len(dates_ERT))))[:len(dates_ERT)]

    # markerline_irr, stemline_irr, baseline_irr = ax.stem(dates_irr, levels_irr,
    #                                          linefmt="--", basefmt="k-")

    # c = ['orange']*len(dates_irr)

    irr_log['color_plot'] = 'orange'
    try:
        for idl in id_left:
            irr_log['color_plot'][idl] = 'green'
    except:
        pass

    ax.bar(dates_irr, levels_irr, color=irr_log['color_plot'])
    # markerline_ERT, stemline_ERT, baseline_ERT = ax.stem(dates_ERT, levels_ERT,
    #                                          linefmt="C4-", basefmt="k-")

    ax.scatter(dates_ERT, levels_ERT, marker='v', color='black')

    # plt.setp(markerline, mec="k", mfc="w", zorder=3)

    # Shift the markers to the baseline by replacing the y-data by zeros.
    # markerline.set_ydata(np.zeros(len(dates)))

    # annotate lines
    vert = np.array(['top', 'bottom'])[(levels_ERT > 0).astype(int)]
    for d, l, r, va in zip(dates_ERT, levels_ERT, names_str, vert):
        ax.annotate(r, xy=(d, l), xytext=(-1, np.sign(l)*3),
                    textcoords="offset points", va=va, ha="right", rotation=25)

    # # remove y axis and spines
    # ax.get_yaxis().set_visible(False)
    # for spine in ["left", "top", "right"]:
    #     ax.spines[spine].set_visible(False)
    ax.set_ylabel('Input water (mL)')
    ax.set_xlabel('Date')

    # ax.margins(y=0.1)
    if show:
        plt.show()
    return ax


def select_from_dates(df, start, end):
    start_datetime = pd.to_datetime(start, format='%d/%m/%Y,%H:%M')
    end_datetime = pd.to_datetime(end, format='%d/%m/%Y,%H:%M')
    mask = (df['datetime'] > start_datetime) & (df['datetime'] <= end_datetime)
    filtered_df = df.loc[mask]
    return filtered_df


def select_from_cycles(df, cycles):
    mask = (df['PRD Cycle nb'].isin(cycles))
    filtered_df = df.loc[mask]
    return filtered_df

# %%
# Load ERT data filenames


def add2pkl_ERT():
    pass


def add2df_ERT(kresipy, selected_files_ERT, ERT_log, TL_flag=False):
    ''' Create a dataframe containing datetimes and the respectives ER values after inversion'''
    df_ERT = []
    for i, f in enumerate(selected_files_ERT):
        date_f = ERT_log['datetime'].loc[ERT_log[ERT_log['Name'] == f].index.to_list()[
            0]]
        kresipy[i].meshResults[0].df['Conductivity(mS/m)']
        if len(df_ERT) == 0:
            df_ERT = pd.DataFrame(
                {date_f: kresipy[i].meshResults[0].df['Conductivity(mS/m)']})
        else:
            df_ERT2merge = pd.DataFrame(
                {date_f: kresipy[i].meshResults[0].df['Conductivity(mS/m)']})
            df_ERT = pd.concat([df_ERT, df_ERT2merge], axis=1)
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


def load_ERT_survey_log(csv2read=('/home/ben/Documents/GitHub/BenjMy/' +
                                  'rhizotron_PRD/imaging_ERT_MALM/' +
                                  'PRD_Measurements_log - 2nd_run.csv'),
                        startDate=None,
                        endDate=None,
                        cycles=[None]
                        ):
    survey_log = pd.read_csv(csv2read,
                             decimal=',',
                             header='infer'
                             )
    survey_log = survey_log.dropna(subset=['Instrument', 'Date'])
    survey_log['Date']
    survey_log['Time (CET in hh:mm)']

    survey_log['datetime_tmp'] = survey_log['Date'] + \
        ' ' + survey_log['Time (CET in hh:mm)']
    # survey_log['datetime'] = pd.to_datetime(survey_log['datetime'], format='%d/%-m/%Y %H:%M')
    survey_log['datetime'] = pd.to_datetime(
        survey_log['datetime_tmp'], format='%d/%m/%Y %H:%M')

    id_ERT = survey_log[survey_log['Name'].str.contains("ERT")].index.tolist()
    id_MALM = survey_log[survey_log['Name'].str.contains(
        "MALM")].index.tolist()

    # survey_log['method'] =
    survey_log["method"] = None
    survey_log['method'][id_ERT] = 'ERT'
    survey_log['method'][id_MALM] = 'MALM'

    survey_log["injection"] = 'Stem'
    id_soil_inj = survey_log[survey_log['Name'].str.contains(
        "MALM2")].index.tolist()
    survey_log['injection'][id_soil_inj] = 'Soil'
    survey_log['injection'][id_ERT] = 'Soil'

    if startDate is not None:
        survey_log = select_from_dates(survey_log, startDate, endDate)

    if type(cycles) == int:
        cycles = [cycles]
    if cycles[0] is not None:
        survey_log = select_from_cycles(survey_log, cycles=cycles)

    return survey_log


def add_process_log(df, f):
    df['inverted_bool'] = np.ones(len(df))*False
    df[df['Name'] == f]['inverted_bool'] = True
    return df


def indiv_invert_ERT(imaging_path, f, recip=5):
    k_indiv = proc.invert_ERT(imaging_path,
                              filename=f,
                              recip=recip
                              )
    proc.plot_ERT(k_indiv, vmin=0, vmax=50,
                  attr="Resistivity(ohm.m)",
                  index=0, ext=['png', 'svg'])
    return k_indiv


def process_ERT(imaging_path, selected_files_ERT, ERT_log, reprocessed=False,
                recip=5):
    k_indiv_merged = []
    for i, f in enumerate(selected_files_ERT):
        print(f)
        date_f = ERT_log['datetime'].loc[ERT_log[ERT_log['Name'] == f].index.to_list()[
            0]]
        if os.path.exists(os.path.join(imaging_path, 'inversionERT') + '/' + f + "bkp.resipy"):
            if reprocessed:
                k_indiv = indiv_invert_ERT(imaging_path, f, recip=recip)
                k_indiv_merged.append(k_indiv)
            else:
                print('reload file from *.resipy backup')
                k_loaded = Project(dirname=imaging_path, typ='R3t')
                k_loaded.loadProject(os.path.join(
                    imaging_path, 'inversionERT') + '/' + f + "bkp.resipy")
                k_indiv_merged.append(k_loaded)
                pass
        else:
            k_indiv = indiv_invert_ERT(imaging_path, f, recip=recip)
            k_indiv_merged.append(k_indiv)
    return k_indiv_merged

def add2df_MALM_icsd(sol_icsd_stk, selected_files_MALM, ERT_log):
    ''' Create a dataframe containing datetimes and the respectives Resistances values '''

    # df_MALM = pd.concat(R_obs_stck, axis=1, keys=[s.name for s in R_obs_stck])
  
    df_MALM_icsd = pd.DataFrame(sol_icsd_stk.T, columns=(selected_files_MALM))
    col_datetime = ERT_log[ERT_log['Name'].isin(selected_files_MALM)]['datetime'].to_list()
    col_injection_type = ERT_log[ERT_log['Name'].isin(selected_files_MALM)]['injection'].to_list()
    tuples = list(zip(df_MALM_icsd.columns.to_list(), 
                      col_datetime,
                      col_injection_type,
                    )
                  )
    index = pd.MultiIndex.from_tuples(tuples, names=["filename", "datetime", "injection"])
    df_MALM_icsd.columns = index
    
    return df_MALM_icsd

def add2df_MALM(R_obs_stck, selected_files_MALM, ERT_log):
    ''' Create a dataframe containing datetimes and the respectives Resistances values '''

    df_MALM = pd.concat(R_obs_stck, axis=1, keys=[s.name for s in R_obs_stck])

    # remove 8 returns files as the sequence does not allow to plot it easily
    # -----------------------------------------------------------------------
    df_MALM = df_MALM.loc[:, ~df_MALM.columns.str.contains('returns')]
    df_MALM.dropna(axis=0, inplace=True)

    # newcols = ERT_log[ERT_log['Name'].isin(selected_files_MALM)]['datetime']
    newcols = ERT_log[ERT_log['Name'].isin(df_MALM.columns)]['datetime']

    df_MALM.columns = newcols

    tuple_cols = []
    cond_same_nb_returns = ~ERT_log['Name'].str.contains('8returns')
    for tc in zip(newcols, ERT_log[(ERT_log['method'] == 'MALM') &
                                   (cond_same_nb_returns)]['injection']):
        tuple_cols.append(tc)
    df_MALM.columns = pd.MultiIndex.from_tuples(tuple_cols)

    return df_MALM


def calc_eucl_dist(pi_EA, x):

    # Distance
    # -----------------------------
    dx = pi_EA[0] - x[0]
    dy = pi_EA[2] - x[2]
    eucl_dist_A_to_x = math.sqrt((dx*dx)+(dy*dy))

    return dx, dy, eucl_dist_A_to_x


def calc_angle(dx, dy):

    # Angle between electrodes
    # -----------------------------
    if dx != 0:
        tan_angle_AB = dy/dx
        if dx < 0:
            angle_AB = np.degrees(np.arctan(tan_angle_AB))
        else:
            angle_AB = np.degrees(np.arctan(tan_angle_AB))-180
    else:
        print('EA and EB same x location')
        angle_AB = 0
    return angle_AB


def sol_polar(xyz_sol):

    pi_EA = np.array([0.25, 0, 0])
    for xs in xyz_sol:
        dx, dy, eucl_dist_A_sol = calc_eucl_dist(pi_EA, xs)
        angle_A_sol = calc_angle(dx, dy)

    return eucl_dist_A_sol, angle_A_sol


def anisotropy(k_MALM):
    seqMALM = k_MALM.surveys[0].df[['a', 'b', 'm', 'n']]
    elecsMALM = k_MALM.surveys[0].elec
    k_MALM.sequence

    old = []
    for e in range(len(elecsMALM)):
        old.append('1 ' + str(e+1))
    new = np.arange(0, len(elecsMALM))

    for i, e in enumerate(old):
        # seqMALM =  seqMALM.replace(e,str(new[i]+1))
        seqMALM = seqMALM.replace(e, new[i]+1)

    angle_AB = []
    angle_A_MN_bary = []

    bary_MN = []
    pi_EA = np.array([0.25, 0, 0])

    eucl_dist_AB = []
    eucl_dist_A_MN_bary = []

    for quad_nb in range(len(seqMALM)):

        # find position of quadripoles
        # -----------------------------
        pi_EB = elecsMALM[['x', 'y', 'z']].iloc[int(
            seqMALM.iloc[quad_nb]['b']-1)].to_numpy()
        pi_EN = elecsMALM[['x', 'y', 'z']].iloc[int(
            seqMALM.iloc[quad_nb]['n']-1)].to_numpy()
        pi_EM = elecsMALM[['x', 'y', 'z']].iloc[int(
            seqMALM.iloc[quad_nb]['m']-1)].to_numpy()

        # Distance between electrodes
        # -----------------------------
        dx = pi_EA[0] - pi_EB[0]
        dy = pi_EA[2] - pi_EB[2]
        eucl_dist_AB.append(math.sqrt((dx*dx)+(dy*dy)))

        bary_MN.append((pi_EM+pi_EN)/2)
        dx_A_MN_bary = pi_EA[0] - bary_MN[quad_nb][0]
        dy_A_MN_bary = pi_EA[2] - bary_MN[quad_nb][2]
        eucl_dist_A_MN_bary.append(
            math.sqrt((dx_A_MN_bary*dx_A_MN_bary)+(dy_A_MN_bary*dy_A_MN_bary)))

        # Angle between electrodes
        # -----------------------------
        if dx != 0:
            tan_angle_AB = dy/dx
            if dx < 0:
                angle_AB.append(np.degrees(np.arctan(tan_angle_AB)))
            else:
                angle_AB.append(np.degrees(np.arctan(tan_angle_AB))-180)

        else:
            print('EA and EB same x location')
            angle_AB.append(0)

        if dx_A_MN_bary != 0:
            tan_angle_A_bary_MN = dy_A_MN_bary/dx_A_MN_bary
            if dx_A_MN_bary < 0:
                angle_A_MN_bary.append(np.degrees(
                    np.arctan(tan_angle_A_bary_MN)))
            else:
                angle_A_MN_bary.append(np.degrees(
                    np.arctan(tan_angle_A_bary_MN))-180)

        else:
            angle_A_MN_bary.append(0)

    bary_MN = np.vstack(bary_MN)
    angle_AB_rad = np.deg2rad(angle_AB)
    angle_A_MN_bary_rad = np.deg2rad(angle_A_MN_bary)

    # Build Dataframe
    # -----------------------------
    df_anisotropy = seqMALM
    df_anisotropy['bary_MNx'] = bary_MN[:, 0]
    df_anisotropy['bary_MNy'] = bary_MN[:, 1]
    df_anisotropy['bary_MNz'] = bary_MN[:, 2]
    df_anisotropy['angle_AB'] = angle_AB
    df_anisotropy['angle_AB_rad'] = angle_AB_rad
    df_anisotropy['eucl_dist_AB'] = eucl_dist_AB
    df_anisotropy['angle_A_MN_bary'] = angle_A_MN_bary
    df_anisotropy['angle_A_MN_bary_rad'] = angle_A_MN_bary_rad
    df_anisotropy['eucl_dist_A_MN_bary'] = eucl_dist_A_MN_bary
    df_anisotropy['resist'] = k_MALM.surveys[0].df['resist']

    # c = plt.scatter(angle_meanR.index, angle_meanR['eucl_dist_AB'].values, c=angle_meanR['resist'].values, s=50, cmap=cm, alpha=0.75,
    #                )

    return df_anisotropy


def normalised_ratio(k_MALM, df_anisotropy):

    obs = k_MALM.surveys[0].df

    def simu_hom(k_MALM):
        # k_MALM_simu_hom = Project(dirname=os.path.join(path + 'inversionMALM', filename), typ='R3t')
        # k_MALM.loadProject(os.path.join(path + 'inversionMALM', filename  + "backup.resipy"))

        k_MALM_simu_hom = k_MALM
        k_MALM_simu_hom.mesh.df['res0'] = 100
        k_MALM_simu_hom.forward()
        # k_MALM_simu_hom.sequence

        return k_MALM_simu_hom

    k_MALM_simu_hom = simu_hom(k_MALM)
    norm_obs = obs['resist']/k_MALM_simu_hom.surveys[0].df['resist']
    df_anisotropy['norm_obs'] = norm_obs

    return norm_obs


def plot_ABMN_pos(k_MALM, seqMALM):

    Pi_EA = np.array([0.25, 0, 0])

    elecsMALM = k_MALM.surveys[0].elec

    Pi_EB = elecsMALM[['x', 'y', 'z']].iloc[seqMALM['b'].unique()-1].to_numpy()
    # Pi_EN = elecsMALM[['x','y','z']].iloc[seqMALM['m'].unique()-1].to_numpy()
    # Pi_EM = elecsMALM[['x','y','z']].iloc[seqMALM['n'].unique()-1].to_numpy()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = plt.scatter(Pi_EA[0], Pi_EA[2], c='red', s=500, alpha=0.75,
                    )
    c = plt.scatter(Pi_EB[:, 0], Pi_EB[:, 2], c='green', s=500, alpha=0.75,
                    )

def plot_anisotropy_polar(df_anisotropy, eucl_dist_A_sol, angle_A_sol, path,
                          norm=True):

    cm = plt.cm.get_cmap('RdYlBu_r')

    key_res = 'resist'
    if norm:
        key_res = 'norm_obs'
    # %% -----------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    # ax.set_thetamin(0)
    # ax.set_thetamax(180)
    # ---- mod here ---- #
    ax.set_theta_zero_location("E")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise

    ax.scatter(np.deg2rad(angle_A_sol), eucl_dist_A_sol,
               c='black', s=50, cmap=cm, alpha=0.75,
               )

    nb_of_angles = len(df_anisotropy['angle_AB_rad'].unique())
    nb_of_dist = len(df_anisotropy['eucl_dist_AB'].unique())

    df_anisotropy_ang_np = []
    for AB_rad in df_anisotropy['angle_AB_rad'].unique():
        cond1 = (df_anisotropy['angle_A_MN_bary_rad'] < AB_rad + 0.02)
        cond2 = (df_anisotropy['angle_A_MN_bary_rad'] > AB_rad - 0.02)
        idvalid = np.where(cond1 & cond2)
        df_anisotropy_ang = df_anisotropy[key_res].iloc[idvalid]
        df_anisotropy_rd = abs(df_anisotropy['eucl_dist_AB'].iloc[idvalid])
        df_anisotropy_ang_np.append([AB_rad*np.ones(len(df_anisotropy_ang)),
                                     df_anisotropy_rd.to_numpy(),
                                     df_anisotropy_ang.to_numpy()])
    df_anisotropy_ang_np = np.hstack(df_anisotropy_ang_np)

    # c = plt.scatter(df_anisotropy_ang_np[0,:], df_anisotropy_ang_np[1,:],
    #                 c=df_anisotropy_ang_np[2,:], s=5, cmap=cm, alpha=0.75,
    #                 )
    len(df_anisotropy_ang_np[0, :])
    c = plt.scatter(df_anisotropy_ang_np[0, :], df_anisotropy_ang_np[1, :],
                    c=df_anisotropy_ang_np[2, :], s=5, cmap=cm, alpha=0.75,
                    vmin=np.percentile(df_anisotropy_ang_np[2, :], 2),
                    vmax=np.percentile(df_anisotropy_ang_np[2, :], 98))

    plt.title('Angle(AB_rad)/dist(AB) = f(resist)')

    # %% -----------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    # ---- mod here ---- #
    ax.set_theta_zero_location("E")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise

    ax.scatter(np.deg2rad(angle_A_sol), eucl_dist_A_sol,
               c='black', s=50, cmap=cm, alpha=0.75,
               )

    nb_of_angles = len(df_anisotropy['angle_AB_rad'].unique())
    nb_of_dist = len(df_anisotropy['eucl_dist_AB'].unique())

    c = plt.scatter(df_anisotropy['angle_AB_rad'], df_anisotropy['eucl_dist_AB'],
                    c=df_anisotropy[key_res], s=5, cmap=cm, alpha=0.75,
                    vmin=np.percentile(df_anisotropy[key_res], 2),
                    vmax=np.percentile(df_anisotropy[key_res], 98))

    plt.title('Angle(AB_rad)/dist(AB) = f(resist)')

    # %% -----------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    # ---- mod here ---- #
    ax.set_theta_zero_location("E")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise

    ax.scatter(np.deg2rad(angle_A_sol), eucl_dist_A_sol,
               c='black', s=50, cmap=cm, alpha=0.75,
               )

    nb_of_angles = len(df_anisotropy['angle_A_MN_bary_rad'].unique())
    nb_of_dist = len(df_anisotropy['eucl_dist_A_MN_bary'].unique())

    c = plt.scatter(df_anisotropy['angle_A_MN_bary_rad'], df_anisotropy['eucl_dist_A_MN_bary'],
                    c=abs(df_anisotropy[key_res]), s=5, cmap=cm, alpha=0.75,
                    vmin=abs(np.percentile(df_anisotropy[key_res], 2)),
                    vmax=abs(np.percentile(df_anisotropy[key_res], 98))
                    )

    # c = plt.scatter(df_anisotropy['angle_A_MN_bary_rad'], df_anisotropy['eucl_dist_A_MN_bary'],
    #                 c=abs(df_anisotropy['resist']), s=5, cmap=cm, alpha=0.75,
    #                 )
    plt.colorbar()

    plt.title('Angle(A-MN)/dist(A-MN) = f(resist)')

    # %% -----------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    # ---- mod here ---- #
    ax.set_theta_zero_location("E")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise

    # angleAB_R_mean = df_anisotropy.groupby(['angle_AB_rad','eucl_dist_A_MN_bary'])[['resist']].mean()
    angleAB_R_mean = df_anisotropy.groupby(['angle_AB_rad'])[
        [key_res, 'eucl_dist_A_MN_bary']].mean()
    distAB_R_mean = df_anisotropy.groupby(['angle_AB_rad'])[[key_res]].mean()

    ax.scatter(np.deg2rad(angle_A_sol), eucl_dist_A_sol,
               c='black', s=50, cmap=cm, alpha=0.75,
               )

    c = plt.scatter(angleAB_R_mean.index, angleAB_R_mean['eucl_dist_A_MN_bary'],
                    c=angleAB_R_mean[key_res], s=5, cmap=cm, alpha=0.75,
                    )
    # c = plt.scatter(angleAB_R_mean.index, angleAB_R_mean['eucl_dist_A_MN_bary'],
    #                 c=angleAB_R_mean['resist'], s=5, cmap=cm, alpha=0.75,
    #                 vmin=np.percentile(angleAB_R_mean['resist'], 2),
    #                 vmax=np.percentile(angleAB_R_mean['resist'], 98)
    #                 )

    plt.title('angleDistAB_f_resist_mean')
    plt.savefig('angleDistAB_f_resist_mean.png', dpi=400)

    # %% -----------------------------------------

    # angleAB_R_mean = df_anisotropy.groupby(['angle_AB_rad','eucl_dist_A_MN_bary'])[['resist']].mean()
    angleAB_R_mean = df_anisotropy.groupby(
        ['angle_AB_rad', 'eucl_dist_A_MN_bary'])[[key_res]].mean()

    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    # ---- mod here ---- #
    ax.set_theta_zero_location("E")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise

    ax.scatter(np.deg2rad(angle_A_sol), eucl_dist_A_sol,
               c='black', s=50, cmap=cm, alpha=0.75,
               )

    for d, new_df in angleAB_R_mean.groupby(level=0):

        anglei = new_df.index.get_level_values(0)
        disti = new_df.index.get_level_values(1)

        # print(anglei[0])
        c = plt.scatter(anglei, disti,
                        c=new_df[key_res], s=5, cmap=cm, alpha=0.75,
                        vmin=np.percentile(df_anisotropy[key_res], 2),
                        vmax=np.percentile(df_anisotropy[key_res], 98)
                        )

    plt.title('angleiDistABi_f_resist')
    plt.savefig('angleiDistABi_f_resist.png', dpi=400)

    return


def petro_plots(scalesurvey,k_indiv_merged,df_SWC,irr_log, ax=None):
    
    if ax == None:
        fig, ax = plt.subplots(figsize=(10, 3))
        
    plot_PRD_effect_SWC(k_indiv_merged,df_SWC,irr_log=irr_log,
                                 unit='m3/m3')
    ax.set_title('SWC variations')
    plt.savefig('../figures/SWC_variations.png', dpi=400)
    
    
    
    fig, ax = plt.subplots(1,sharex=True)
    
    color = 'tab:red'
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y1-axis', color = color)
    ax2 = ax.twinx()
    
    color = 'tab:green'
    ax2.set_ylabel('Y2-axis', color = color)
            
    # fig, axs = plt.subplots(3,sharex=False)
    
    # ax_irr = surveyPRD.plot_timeline(ERT_log,irr_log,ax=axs[0],dropMALM=True)
    ax_scale = plot_weight_dynamic(scalesurvey,ax=ax, color='red')
    
    ax_SWC = plot_PRD_effect_SWC(k_indiv_merged,df_SWC,irr_log=irr_log,ax=ax2)
    ax.set_title('ER VS weight data')
    plt.savefig('../figures/SWC_weight_LR.png', dpi=400)
    
    #%%
    
    fig, ax = plt.subplots(1,sharex=True)
    color = 'tab:red'
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y1-axis', color = color)
    ax2 = ax.twinx()
    color = 'tab:green'
    ax2.set_ylabel('Y2-axis', color = color)
            
    ax_scale = plot_weight_dynamic(scalesurvey,ax=ax, color='red')
    
    ax_SWC = plot_PRD_effect_SWC(k_indiv_merged,df_SWC,irr_log=irr_log,ax=ax2,
                                           LR=False, topdown=True)
    ax.set_title('ER VS weight data')
    plt.savefig('../figures/SWC_weight_topdown.png', dpi=400)
    
    
    #%%
    
    fig, ax = plt.subplots(1,sharex=True)
    color = 'tab:red'
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y1-axis', color = color)
    ax2 = ax.twinx()
    color = 'tab:green'
    ax2.set_ylabel('Y2-axis', color = color)
            
    ax_scale = plot_weight_dynamic(scalesurvey,ax=ax, color='red')
    
    ax_SWC = plot_PRD_effect_SWC(k_indiv_merged,df_SWC,irr_log=irr_log,ax=ax2,
                                           LR=False, topdown=True)
    ax.set_title('ER VS weight data')
    plt.savefig('../figures/SWC_weight_mean.png', dpi=400)
    
    pass

# %%


def load_irr_log_drive(startDate=None,
                       endDate=None,
                       cycles=[None]):
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
    irr_log['datetime'] = pd.to_datetime(
        irr_log['datetime_tmp'], format='%d/%m/%Y %H:%M')

    if startDate is not None:
        irr_log = select_from_dates(irr_log, startDate, endDate)

    if type(cycles) == int:
        cycles = [cycles]
    if cycles[0] is not None:
        cycles.insert(0, min(cycles)-1)
        irr_log = select_from_cycles(irr_log, cycles=cycles)

    return irr_log

# Load scale data filenames


def plot_weight_dynamic(scalesurvey, ax=None, **kwargs):

    if ax == None:
        fig, ax = plt.subplots(figsize=(18, 5))
       
    # def date2yday(x):
    # """Convert matplotlib datenum to days since 2018-01-01."""
    # y = x - mdates.date2num(datetime.datetime(2018, 1, 1))
    # return y

    color = 'black'
    if 'color' in kwargs:
        color = kwargs['color']
        
    # def yday2date(x):
    #     """Return a matplotlib datenum for *x* days after 2018-01-01."""
    #     y = x + mdates.date2num(datetime.datetime(2018, 1, 1))
    #     return y
    
    
    # secax_x = ax.secondary_xaxis('top', functions=(date2yday, yday2date))
    # secax_x.set_xlabel('yday [2018]')

    scalesurvey.plot('datetime', 'weight (kg)', style='.', ax=ax,
                     color=color)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    # ax.xaxis.set_minor_locator(mdates.HourLocator(interval=4))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
    # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))

    # Set title and labels for axes
    ax.set(xlabel="Date",
           ylabel="weight (kg)",
           title="Raw data sheet")
    
    return ax


def load_scale_data(startDate=None,
                    endDate=None,
                    cycles=[None]):

    sheet_id = '1GqDKU8fC0Pw_DQpUD-MBLpN-AXWrAzLuO1hABxxYIhc'
    sheet_name = ['apr_6-11', 'apr_11-13', 'apr_13-21', 'apr_21-26', 'apr_28-may_2',
                  'may_2-9', 'may_9-12', 'may_13-16', 'may_16-19', 'may_19-23',
                  'may_23-31','may_31-june8','jun_8-15','jun_15-22','jun_23-29',
                  'jun_29-jul_5','jul_5-11',
                  ]

    root_name = 'transpiration_track_scale_II.xlsx - '
    
    list_data_scale = []
    for sn in sheet_name:
        # url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sn}'
        data_scale_i = pd.read_csv('../load_cells/rawData/' + root_name + sn +'.csv', decimal=',')

        data_scale_i['abs_date'] = pd.to_datetime(data_scale_i['date'] + " " + data_scale_i['time'],
                                                  infer_datetime_format=True)
        datestr = data_scale_i['date'] + " " + data_scale_i['time']
        data_scale_i['abs_date'] = pd.to_datetime(data_scale_i['date'] + " " + data_scale_i['time'],
                                                  format="%d/%m/%Y %H:%M")
        initial_date = data_scale_i['abs_date'][0]
        dates = []
        for i in range(len(data_scale_i['time'])):
            dates.append(initial_date +
                         pd.Timedelta(seconds=data_scale_i['sec'][i]))
        date_time = [t.strftime("%Y-%m-%d %H:%M:%S") for t in dates]
        data_scale_i['datetime'] = date_time
        data_scale_i['datetime'] = pd.to_datetime(data_scale_i['datetime'])
        list_data_scale.append(data_scale_i)

    df = pd.concat(list_data_scale, axis=0, ignore_index=False)

    if startDate is not None:
        df = select_from_dates(df, startDate, endDate)
    # if cycles[0] is not None:
    #     cycles.insert(0, min(cycles)-1)
    #     df = select_from_cycles(df, cycles=cycles)

    return df
            
        
        
def Archie_rho2sat(rho, rFluid, porosity, a=1.0, m=2.0, n=2.0):
    # See petro file pyCATHY  wrapper envs for chanages
    '''
    rho: resistivity
    𝑆𝑤 : water saturation
    𝜙: the porosity of the soil
    𝜎_{𝑤} is the conductivity of the pore fluid
    𝑎, 𝑚, and 𝑛 are empirically derived parameters
    𝑎 is the tortuosity factor
    𝑚 is the cementation exponent
    𝑛 is the saturation exponent
    Returns
    -------
    𝑆𝑤 : water saturation
    '''
    return (rho/rFluid/porosity**(-m))**((-1/n))

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
