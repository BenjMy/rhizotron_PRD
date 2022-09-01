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
import surveyPRD

imaging_path = '../imaging_ERT_MALM/'

# %%


def get_cmd():
    parse = argparse.ArgumentParser()

    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[None], required=False) #[0,1,2,3,4,5]
    period.add_argument('-cycle', '--cycle', nargs='+',
                        help='list of cycle', type=int, default=[3,4,5,6,7,8,9], required=False)
    # period.add_argument('-cycle', '--cycle', nargs='+',
    #                     help='list of cycle', type=int, default=[4,5,6,7,8,9], required=False)
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
                               help='Apply Archie trans.', default=0, required=False)
    process_param.add_argument('-pareto', type=int,
                               help='pareto', default=0, required=False)
    process_param.add_argument('-wr', type=float,
                               help='reg. weight', default=1, required=False)
    args = parse.parse_args()
    return(args)


# %%
# print(args)
# args.cycle = [1,2,3,4]
ERT_log_all = surveyPRD.load_ERT_survey_log()
irr_log_all = surveyPRD.load_irr_log_drive()
ax = surveyPRD.plot_timeline(ERT_log_all, irr_log_all, dropMALM=True)
# tdelta = (ERT_log_all['datetime'].iloc[0]-ERT_log_all['datetime'].iloc[-1])/10
# ax.get_xaxis().set_major_locator(mdates.DayLocator(interval=5))
# ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))


# %%
# ERT_log['PRD Cycle nb']
# args.cycle

# %%
# k.createBatchSurvey(imaging_path + 'filenames_ERT')
# k.invert(parallel=True)
def run_ERT_ind(args, selected_files_ERT, ERT_log):
    k_indiv_merged = surveyPRD.process_ERT(imaging_path,
                                           selected_files_ERT,
                                           ERT_log,
                                           reprocessed=bool(args.reprocessed),
                                           recip=args.recErr,
                                           )
    
    for fi in range(len(selected_files_ERT)):
        proc.plot_ERT(k_indiv_merged[fi], vmin=0, vmax=2,
                      attr="Resistivity(log10)")
        proc.plot_ERT(k_indiv_merged[fi], vmin=0, vmax=50,
                      attr="Resistivity(ohm.m)")
    
    # for i, ki in enumerate(k_indiv_merged):
        # nb_of_rejected_quad[i] = ki.surveys[0]
        # final_rms = 
        # nb_of_iter = 
    
    # %%
    # surveyPRD.add2pkl_ERT()
    df_ERT = surveyPRD.add2df_ERT(k_indiv_merged,
                                  selected_files_ERT,
                                  ERT_log)
    
    # surveyPRD.df_ERT_diff(df_ERT,selected_files_ERT,
    #                 background_diff=False)
    
    # def df_ERT_diff(df_ERT,selected_files_ERT,
    #                 background_diff=False):
    
    #     df_ERT['diff_level'] = False
    #     df_ERT = df_ERT.set_index('diff_level', append=True).unstack('diff_level')

           
    #     # df_diff = pd.DataFrame()
    #     for i, f in enumerate(selected_files_ERT):
    #         if (i < len(selected_files_ERT)-1):
    #             fig, ax = plt.subplots(
    #                                    constrained_layout=False)
                
    #             d1=ERT_log[ERT_log['Name'] == selected_files_ERT[i+1]]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
    #             if background_diff:
    #                 id_file = 0
    #                 d0=ERT_log[ERT_log['Name'] == selected_files_ERT[0]]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
    #             else:
    #                 id_file = i
    #                 d0=ERT_log[ERT_log['Name'] == selected_files_ERT[0]]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
    
        
        
        
    # %%
    fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))
    ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)
    ax = surveyPRD.plot_PRD_effect_ER(
        k_indiv_merged, df_ERT, irr_log, ax=ax[1])
    plt.savefig(imaging_path+'inversionERT/' + 'PRDeffect' + str([*args.cycle]),
                dpi=450)
    
    if args.petro:
        rFluid = 1
        porosity = 0.4
        df_SWC = df_ERT.copy()
        for cols in list(df_ERT.columns):
            if type(cols)==str:
                pass
            else:
                df_SWC[cols] = surveyPRD.Archie_rho2sat(df_ERT[cols].to_numpy(),
                                                          rFluid, porosity, 
                                                          a=1.0, m=2.0, n=2.0)
            
        
    return k_indiv_merged, df_ERT


def run_TL(selected_files_ERT):
    idTL = []
    for sf in selected_files_ERT:
        idTL.append(ERT_log[ERT_log['Name'] == sf]['PRD Cycle nb'].tolist()[0])

    # selected_files_ERT_reverse = list(reversed(selected_files_ERT))
    # idTL_reverse = list(reversed(idTL))
    cycle_ERT_time = list(np.unique(idTL))

    # regType=1 # background constrainst inv

    k_TL = proc.invert_ERT_TL(
        imaging_path,
        files=selected_files_ERT,
        regType=1,
        recip=5,
        idfileNames=cycle_ERT_time,
        reprocessed=bool(args.reprocessed),
        # reprocessed=True,
    )
    
    for fi in range(len(selected_files_ERT)):
        # proc.plot_ERT(k_TL[0], vmin=-10, vmax=1,
        #               attr="Sensitivity_map(log10)", index=fi)
        proc.plot_ERT(k_TL[0], vmin=0, vmax=2,
                      attr="Resistivity(log10)", index=fi)
        proc.plot_ERT(k_TL[0], vmin=0, vmax=50,
                      attr="Resistivity(ohm.m)", index=fi)
        if fi > 0:
            proc.plot_ERT(k_TL[0], vmin=-20, vmax=20,
                          attr="difference(percent)", index=fi)
    
    
    # k_TL = proc.invert_ERT_TL(
    #     imaging_path,
    #     files=selected_files_ERT,
    #     regType=2,
    #     recip=5,
    #     idfileNames=cycle_ERT_time,
    #     reprocessed=bool(args.reprocessed),
    #     # reprocessed=True,
    # )
    
    # for fi in range(len(selected_files_ERT)):
    #     proc.plot_ERT(k_TL[0], vmin=-10, vmax=1,
    #                   attr="Sensitivity_map(log10)", index=fi)
    #     proc.plot_ERT(k_TL[0], vmin=0, vmax=2,
    #                   attr="Resistivity(log10)", index=fi)
    #     proc.plot_ERT(k_TL[0], vmin=0, vmax=50,
    #                   attr="Resistivity(ohm.m)", index=fi)
    #     if fi > 0:
    #         proc.plot_ERT(k_TL[0], vmin=-20, vmax=20,
    #                       attr="difference(percent)", index=fi)
    
    

    # %%
    # len(k_TL[0].surveys)   
    # pl = pv.Plotter(shape=(1,len(k_TL[0].surveys)-1), window_size=[800, 400],
    #                 notebook=True)
    # # pl.show_bounds(font_size=10)
    
    # for i in range(len(k_TL[0].surveys)-1):
    #     print(i)
    #     pl.subplot(0,i)
    #     mesh =  k_TL[0].mesh.df
    #     if i > 0:
    #         k_TL[0].showResults(
    #                       index=i,
    #                       attr="difference(percent)",
    #                       ax=pl, 
    #                       vmin=-20, vmax=20,
    #                       color_map='jet',
    #                       background_color='white',
    #                       pvgrid = True,
    #                       use_pyvista=True,
    #                       pvshow=False,
    #                       xlim=[mesh['X'].min(),mesh['X'].max()],
    #                       ylim=[mesh['Y'].min(),mesh['Y'].max()], 
    #                       zlim=[mesh['Z'].min(),mesh['Z'].max()], 
    #                       )
    #         # pl.x_axis.tick_size += 10
    # pl.screenshot(imaging_path + 'foo.png')
    # pl.show() 

        # pl.save_graphic("img.svg")



# %%


def run_icsd(selected_files_MALM, selected_files_ERT, k_indiv_merged,
             inversionPathMALM):
    ''' Plot the observed data '''
    
    k_MALM = []
    f_MALM = []
    R_obs_stck = []
    fig = plt.figure()
    ax = fig.add_subplot()
       
    for i, f in enumerate(selected_files_MALM):
        background_ERT_time = int(f.split('_')[2])
        for i, n in enumerate(ERT_log['Name'][ERT_log['method'] == 'ERT']):
            if 'PRD_ERT_' + str(background_ERT_time) in n:
                index_ERT_backgrd = i

        R_obs = proc.prepare_MALM(imaging_path+inversionPathMALM,
                                  f,
                                  k_indiv_merged[index_ERT_backgrd],
                                  filter_seq=True,
                                  filter_seq_rec=False,
                                  percent_rec=1e9,
                                  ax=ax
                                  )
        if len(R_obs_stck)>0:
            if len(R_obs)==len(R_obs_stck[0]):
                R_obs_stck.append(R_obs)
                f_MALM.append(f)
        else:
            R_obs_stck.append(R_obs)
            f_MALM.append(f)
    
    for ff in f_MALM:
        if'return' in ff:
            f_MALM.remove(ff)


    #%%
    df_MALM = surveyPRD.add2df_MALM(R_obs_stck,
                                    f_MALM,
                                    ERT_log)
    
    #%%
    fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))
    
    # df_SWC
    
    ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)

    ax = surveyPRD.plot_PRD_effect_MALM(k_indiv_merged, df_MALM, ax=ax[1])
    # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))

    plt.xticks(rotation=35)
    plt.tight_layout()

    if type(args.cycle) is int:
        args.cycle = [args.cycle]
    plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect' + str([*args.cycle]),
                dpi=450)
    plt.close('all')


    #%%

    proc.plot_scatter_MALM(imaging_path+inversionPathMALM,
                           df_MALM,
                           k_indiv_merged,
                           ERT_log,
                           vmin=-125, vmax=550,
                           )
    plt.close('all')


    #%%
    df_MALM_diff_b = proc.plot_scatter_MALM_diff(imaging_path+inversionPathMALM,
                                df_MALM,
                                k_indiv_merged,
                                f_MALM,
                                ERT_log,
                                background_diff=True,
                                vmin=25, vmax=150,
                                )

    fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))   
    ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)
    ax = surveyPRD.plot_PRD_effect_MALM_diff(k_indiv_merged, df_MALM_diff_b, ax=ax[1])
    # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))
    plt.xticks(rotation=35)
    plt.tight_layout()

    if type(args.cycle) is int:
        args.cycle = [args.cycle]
    plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect_diff_background' + str([*args.cycle]),
                dpi=450)
    plt.close('all')

    
    
    df_MALM_diff = proc.plot_scatter_MALM_diff(imaging_path+inversionPathMALM,
                                df_MALM,
                                k_indiv_merged,
                                f_MALM,
                                ERT_log,
                                background_diff=True,
                                vmin=25, vmax=150,
                                )
    
    fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))   
    ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)
    ax = surveyPRD.plot_PRD_effect_MALM_diff(k_indiv_merged, df_MALM_diff, ax=ax[1])
    # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))
    plt.xticks(rotation=35)
    plt.tight_layout()

    if type(args.cycle) is int:
        args.cycle = [args.cycle]
    plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect_diff' + str([*args.cycle]),
                dpi=450)
    plt.close('all')

    
    
    #%%
    if len(f_MALM) == 2*len(k_indiv_merged):
        df_MALM_diff_StemSoil = proc.plot_scatter_MALM_diff_stem_soil(
                                            imaging_path+inversionPathMALM,
                                            df_MALM,
                                            k_indiv_merged,
                                            f_MALM,
                                            ERT_log,
                                            vmin=25, vmax=1e4,
                                            )
    plt.close('all')

        
        

    #%%
    ERT_log
    k_MALM = []
    for i, f in enumerate(selected_files_MALM):
        print('*'*24)
        print(f)
        print('*'*24)

        background_ERT_time = int(f.split('_')[2])
        for j, n in enumerate(ERT_log['Name'][ERT_log['method'] == 'ERT']):
            if 'PRD_ERT_' + str(background_ERT_time) in n:
                index_ERT_backgrd = j
        # ERT_log[ERT_log['Name'].isin(['PRD_ERT_' + str(background_ERT_time)])]
        # ERT_log[ERT_log['Name']==('PRD_ERT_' + str(background_ERT_time))]
        # print(background_ERT_time)
        outMALM = proc.prepare_icsd(imaging_path + inversionPathMALM,
                                        f,
                                        k_indiv_merged[index_ERT_backgrd],
                                        reprocessed=bool(args.reprocessed),
                                        # reprocessed=True,
                                        nVRTe_rows=9, nVRTe_cols=9,
                                        filter_seq_rec=args.filter_seq_rec,
                                        filter_seq=args.filter_seq,
                                        reduce2d=True,
                                        )
        k_MALM.append(outMALM[0])
        
    [_, imin, R_obs, nodes, R_icsd] = outMALM
    plt.close('all')


    # %%

    m0 = []
    # selected_files_MALM_without_8_returns
    j = 0
    for i, f in enumerate(selected_files_MALM):
        if '8returns' in f:
            continue
        # j = ind_MALM_SOIL_STEM(i,selected_files_MALM,k_indiv_merged)
        m0i,ax,fig = proc.m0_MALM(imaging_path + inversionPathMALM + f,
                               method_m0='F1', typ=args.dim,
                               show=True
                               )
        fig.savefig(os.path.join(imaging_path,inversionPathMALM, f,'m0.png'), dpi=300)
        m0.append(m0i)
        pl = proc.plot_m0_MALM(k_MALM[j], m0[j], pl=None, ext=['png'], show=False,
                          index=0)
        pl.close()
        pv.close_all()


        j += 1
        plt.close('all')

       
        # m0.append(proc.m0_MALM(imaging_path + inversionPathMALM + f,
        #                        method_m0='Pearson', typ='2d')
        #           )
        # proc.plot_m0_MALM(k_MALM[i], m0[i], pl=None, ext=['png'], show=False,
        #                   index=0)
        
        # plt.close('all')

    # %%

    sol_stk = []
    f_names = []
    j = 0

    for i, f in enumerate(selected_files_MALM):
        if (('8returns' in f) and (args.filter_seq==True)):
            continue
        
        print('*'*24)
        print(f)
        print('*'*24)
        pareto = bool(args.pareto)
        prior = False
        # for prior in [True,False]:
        sol, ax, fig = proc.invert_MALM(imaging_path + inversionPathMALM + f,
                                    wr=args.wr, typ=args.dim, pareto=pareto, show=True,
                                    prior = prior,
                                    # fname_sim='VRTeSim_icsd.txt'
                                    )
        sol_stk.append(sol)
        fig.savefig(os.path.join(imaging_path,inversionPathMALM, f,'icsd.png'), dpi=300)
        plt.close('all')
        # proc.plot_MALM(k_MALM[i], nodes, imin, sol[i])
        if pareto:
            sol = sol_stk[j][0].solution
            pl = proc.plot_MALM(k_MALM[j], sol, show=False, ext=['png'],
                            prior = prior)
            pl.close()
            # proc.plot_MALM(k_MALM[i],sol,show=True,ext=['svg'])
        else:
            pl = proc.plot_MALM(k_MALM[j], sol_stk[j], show=False, ext=['png'],
                            prior = prior)
            # proc.plot_MALM(k_MALM[i],sol,show=True,ext=['svg'])
            pl.close()
        j += 1
        f_names.append(f)
        pv.close_all()




    if pareto:
        sol_icsd_stk_new = [x[0].solution.x for x in sol_stk]
        sol_icsd_stk_new = np.vstack(sol_icsd_stk_new)
    else:
        sol_icsd_stk_new = [x.x for x in sol_stk]
        sol_icsd_stk_new = np.vstack(sol_icsd_stk_new)

    df_MALM_icsd = surveyPRD.add2df_MALM_icsd(sol_icsd_stk_new,
                                              f_names,
                                              ERT_log)
    
    
    df_MALM_icsd.to_csv(imaging_path + 
                        inversionPathMALM + 
                        'df_MALM_icsd' + 
                        str([*args.cycle]) 
                        + '.csv', index=False)
    
    np.savetxt(imaging_path + 
                inversionPathMALM + 
                'vrte_nodes_in_mesh.txt',
                imin)
    
    
    #%%

    
    fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))   
    ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)
    
    df_MALM_icsd = surveyPRD.plot_PRD_effect_icsd(k_indiv_merged, 
                                                  imin,
                                                  df_MALM_icsd, irr_log,
                                                  ax=ax[1],
                                                  hours_interval=10)
    # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))
    plt.xticks(rotation=35)
    plt.tight_layout()

    if type(args.cycle) is int:
        args.cycle = [args.cycle]
    plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect_icsd' + str([*args.cycle]),
                dpi=450)
    
    
# %%
if __name__ == '__main__':

    '''
    - Read files
    - Invert ERT
    - Invert MALM
    '''

    # Read and parse args
    # -----------------------------------
    args = get_cmd()
    # args.startD = '8/6/2022,9:59'
    # args.endD = '8/6/2022,11:41'    
    # args.startD = '21/6/2022,13:50'
    # args.endD = '26/6/2022,14:50'    
    
    # args.startD = '5/7/2022,13:50'
    # args.endD = '9/7/2022,14:50'
    
    

    inversionPathMALM = surveyPRD.definePaths(args)

    # Read and parse log files
    # -----------------------------------
    f_ERT, f_MALM,  ERT_log, irr_log = surveyPRD.load_log_files(args)

    # ERT inversion
    # -----------------------------------
    k_indiv_ERT, df_ERT = run_ERT_ind(args, f_ERT, ERT_log)

    # TL ERT inversion
    # -----------------------------------
    if args.TL:
        run_TL(f_ERT)

    # %% MALM icsd inversion
    # -----------------------------------
    if args.icsd:
        run_icsd(f_MALM, f_ERT, k_indiv_ERT, inversionPathMALM)
