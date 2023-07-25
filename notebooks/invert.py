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

imaging_path = '../imaging_ERT_MALM_PaperREV1/'

# %%


def get_cmd():
    parse = argparse.ArgumentParser()

    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[None], required=False) #[0,1,2,3,4,5]
    period.add_argument('-cycle', '--cycle', nargs='+',
                        # help='list of cycle', type=int, default=[3,4,5,6,7,8], required=False)
                        # help='list of cycle', type=int, default=[0,1,2,3,4,5,6,7,8], required=False)
                        # help='list of cycle', type=int, default=[3,4,5,6,7,8], required=False)
                        help='list of cycle', type=int, default=[6,7], required=False)
                        # help='list of cycle', type=int, default=[5,6], required=False)
                        # help='list of cycle', type=int, default=[-99], required=False)
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
        '-TL', type=int, help='TimeLapse', default=1, required=False)
    process_param.add_argument(
        '-TLreg', type=int, help='TimeLapse reg mode', default=1, required=False)
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
    process_param.add_argument('-pareto', type=int,
                               help='pareto', default=0, required=False)
    process_param.add_argument('-wr', type=float,
                               help='reg. weight', default=10, required=False)
    
    process_param.add_argument('-anisotropy', type=float,
                               help='anisotropy procesing', default=0, required=False)
    args = parse.parse_args()
    
    
    # args.startD = '21/6/2022,13:50'
    # # args.startD = '29/6/2022,09:29'
    # args.endD = '29/6/2022,10:13'
    
    # args.startD = '11/7/2022,11:19' 
    # args.endD = '12/7/2022,13:05'

    # args.startD = '11/7/2022,11:19' 
    # args.endD = '12/7/2022,13:05'
    
    
    # args.startD = '29/6/2022,9:00'
    # args.endD = '5/7/2022,17:00'
    
    # ALL
    # 13/5/2022	16:25
    # 12/7/2022	12:50

    # args.startD = '29/6/2022,09:29'
    # args.endD = '5/7/2022,16:35'
    
    # LAST
    # args.startD = '11/7/2022,15:49'
    # args.endD = '11/7/2022,17:16'
       
    return(args)

# %%
# print(args)
# args.cycle = [1,2,3,4]
ERT_log_all = surveyPRD.load_ERT_survey_log()
irr_log_all = surveyPRD.load_irr_log()
ax = surveyPRD.plot_timeline(ERT_log_all, irr_log_all, dropMALM=True)
# tdelta = (ERT_log_all['datetime'].iloc[0]-ERT_log_all['datetime'].iloc[-1])/10
# ax.get_xaxis().set_major_locator(mdates.DayLocator(interval=5))
# ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))

print(ERT_log_all)

# %%
# ERT_log['PRD Cycle nb']
# args.cycle

# %%
# k.createBatchSurvey(imaging_path + 'filenames_ERT')
# k.invert(parallel=True)
def run_ERT_ind(
                args, 
                selected_files_ERT, 
                ERT_log,
                **kwargs,
                ):
    k_indiv_merged = surveyPRD.process_ERT(imaging_path,
                                           selected_files_ERT,
                                           ERT_log,
                                           reprocessed=bool(args.reprocessed),
                                           recip=args.recErr,
                                           )
    
    df_rms = []
    for fi in range(len(selected_files_ERT)):
        proc.plot_ERT(k_indiv_merged[fi], vmin=0, vmax=2,
                      attr="Resistivity(log10)",fi=selected_files_ERT[fi])
        proc.plot_ERT(k_indiv_merged[fi], vmin=0, vmax=50,
                      attr="Resistivity(ohm.m)",fi=selected_files_ERT[fi])
    
    # for i, ki in enumerate(k_indiv_merged):
        # nb_of_rejected_quad[i] = ki.surveys[0]
        # final_rms = 
        # nb_of_iter = 

        
    # Plot evolution RMS
    # ---------------------------------------------------------------
    df_rms = (proc.getR3out(imaging_path,
                           selected_files_ERT,
                           )
                  )
    
    df_rms['datetime'] = ''
    for i, dfn in enumerate(df_rms.name):
        test = ERT_log[ERT_log['Name']==dfn]['datetime'].values[0]
        df_rms['datetime'].iloc[i] = test
    
    # df_rms['datetime'] 
    df_rms.set_index('datetime', inplace=True)
    
    
    mosaic = '''A
                B
            '''
    fig, axs = plt.subplot_mosaic(mosaic,sharex=True,figsize=(7,4))
    
    df_rms.plot.hist(y='resRMS',ax=axs['A'])
    df_rms.plot.hist(y='read',ax=axs['B'])
    plt.savefig(imaging_path+'inversionERT/' + 'performance.png', dpi=350)
    
    df_rms.to_csv(imaging_path+'inversionERT/RMS' + str([*args.cycle])+'.csv')

    # %%
    # surveyPRD.add2pkl_ERT()
    df_ERT = surveyPRD.add2df_ERT(k_indiv_merged,
                                  selected_files_ERT,
                                  ERT_log)
      
    
    df_ERT.to_csv(imaging_path+'inversionERT/dfERT' + str([*args.cycle])+'.csv')

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
        k_indiv_merged, df_ERT, irr_log, ax=ax[1],**kwargs)
    # ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%B %d'))
    plt.savefig(imaging_path+'inversionERT/' + 'PRDeffect' + str([*args.cycle]),
                dpi=450)
    
    if args.petro:
                
        from datetime import datetime
        tmin = datetime.strptime('08/06/2022 12:29:00', '%d/%m/%Y %H:%M:%S')
        tmax = datetime.strptime('15/06/2022 16:19:00', '%d/%m/%Y %H:%M:%S')

        porosity = 0.55
        df_SWC = df_ERT.copy()
        for cols in list(df_ERT.columns[:-1]):
            
            sw_conductivity = 2470*(1e-6/1e-2)
            if (tmin < cols) & (tmax > cols):
                sw_conductivity = 2000*(1e-6/1e-2)
                print('Cycle 3 tap water detected - change water conductivity')
            rFluid = 1/sw_conductivity

            if type(cols)==str:
                pass
            else:
                ERi = 1/(df_ERT[cols]*1e-3) # Conductivity(mS/m)
                # print(cols)
                # print(ERi.mean())
                # df_SWC[cols].mean()
                # min(ERi)
                df_SWC[cols] = surveyPRD.Archie_rho2sat(ERi,
                                                        rFluid, porosity, 
                                                        a=1.0, m=1.9, n=1.2)
                df_SWC[cols] = df_SWC[cols]*porosity
        df_SWC.to_csv(imaging_path+'inversionERT/dfSWC' + str([*args.cycle])+'.csv')

        
    return k_indiv_merged, df_ERT


def run_TL(selected_files_ERT):
    idTL = []
    for sf in selected_files_ERT:
        idTL.append(ERT_log[ERT_log['Name'] == sf]['PRD Cycle nb'].tolist()[0])

    # selected_files_ERT_reverse = list(reversed(selected_files_ERT))
    # idTL_reverse = list(reversed(idTL))
    cycle_ERT_time = list(np.unique(idTL))
    print('***'*12)
    print(cycle_ERT_time)
    print('***'*12)

    # regType=1 # background constrainst inv

    k_TL = proc.invert_ERT_TL(
        imaging_path,
        files=selected_files_ERT,
        regType=args.TLreg,
        recip=5,
        idfileNames=cycle_ERT_time,
        reprocessed=bool(args.reprocessed),
        # reprocessed=True,
    )
    
    for fi in range(len(selected_files_ERT)):
        # proc.plot_ERT(k_TL[0], vmin=-10, vmax=1,
        #               attr="Sensitivity_map(log10)", index=fi)
        proc.plot_ERT(k_TL[0], vmin=0, vmax=2,
                      attr="Resistivity(log10)", index=fi,
                      fi=selected_files_ERT[fi])
        proc.plot_ERT(k_TL[0], vmin=0, vmax=50,
                      attr="Resistivity(ohm.m)", index=fi,
                      fi=selected_files_ERT[fi],
                      )
        if fi > 0:
            proc.plot_ERT(k_TL[0], vmin=-20, vmax=20,
                          attr="difference(percent)", index=fi,
                          showElec=False, fi=selected_files_ERT[fi])
    
    
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
             inversionPathMALM,
             **kwargs):
    ''' Plot the observed data '''
    

    k_MALM = []
    f_MALM = []
    R_obs_stck = []
    fig = plt.figure()
    ax = fig.add_subplot()
       
    for i, f in enumerate(selected_files_MALM):
        if '8returns' in f:
            continue
        background_ERT_time = int(f.split('_')[2])
        for i, n in enumerate(ERT_log['Name'][ERT_log['method'] == 'ERT']):
            if 'PRD_ERT_' + str(background_ERT_time) in n:
                index_ERT_backgrd = i

        R_obs = proc.prepare_MALM(imaging_path+inversionPathMALM,
                                  f,
                                  k_indiv_merged[index_ERT_backgrd],
                                  filter_seq=True,
                                  filter_seq_rec=False,
                                  percent_rec=10,
                                  ax=ax,
                                  m=71
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
    
    
    #%% Anisotropy 
    if bool(args.anisotropy):
        k_MALM_anisotropy = []
        m0_anisotropy = []
        for i, f in enumerate(selected_files_MALM):
            if '8returns' in f:
                # b_return = [1,4,8,17,20,24,33,36,40,49,52,64]
                b_return = [1,8,17,24,36,49,64]

                for k in range(len(b_return)):
                    print('A'*24)
                    print(k,f)
                    print('A'*24)
                    fnew = f.split('.csv')[0] + '_B'+ str(b_return[k])
                    background_ERT_time = int(f.split('_')[2])
                    for j, n in enumerate(ERT_log['Name'][ERT_log['method'] == 'ERT']):
                        if 'PRD_ERT_' + str(background_ERT_time) in n:
                            index_ERT_backgrd = j
                    outMALM = proc.prepare_icsd(imaging_path + inversionPathMALM,
                                                    fnew,
                                                    k_indiv_merged[index_ERT_backgrd],
                                                    # reprocessed=1, #bool(args.reprocessed),
                                                    reprocessed=bool(args.reprocessed),
                                                    nVRTe_cols=6*2+1, nVRTe_rows=4*2+1,
                                                    filter_seq_rec=args.filter_seq_rec,
                                                    filter_seq=args.filter_seq,
                                                    b=b_return[k],
                                                    reduce2d=True,
                                                    )
                    k_MALM_anisotropy.append(outMALM[0])
                    
                    m0i,ax,fig = proc.m0_MALM(imaging_path + inversionPathMALM + fnew,
                                           method_m0='F1', typ=args.dim,
                                           show=True,
                                           retElec=k_indiv_merged[0].elec[['x','z']].iloc[b_return[k]-1].to_numpy()
                                           )
                    plt.close()
                    fig.savefig(os.path.join(imaging_path,inversionPathMALM, fnew,'m0' + '_breturn' + str(b_return[k]) + '.png')
                                , dpi=300)
                    m0_anisotropy.append(m0i)
                    
                    
                    icsd, sol, ax, fig = proc.invert_MALM(imaging_path + inversionPathMALM + fnew,
                                                wr=args.wr, typ=args.dim, pareto=args.pareto, show=True,
                                                prior = False,                                                        
                                                clim=[0,0.1],
                                                retElec=k_indiv_merged[0].elec[['x','z']].iloc[b_return[k]-1].to_numpy()
                                                # fname_sim='VRTeSim_icsd.txt'
                                                )
                    


                    # sol_anisotropy_stk.append(sol)
                    # icsd_stk.append(icsd)
                    fig.savefig(os.path.join(imaging_path,inversionPathMALM, fnew,'icsd.png'), dpi=300)
                    plt.close('all')
                    
                    # pl = proc.plot_m0_MALM(k_MALM_anisotropy[k], m0_anisotropy[k], pl=None, ext=['png'], show=False,
                    #                   index=0)
                    
                # [_, imin, R_obs, nodes, R_icsd] = outMALM
                plt.close('all')
            
                # plt.plot(R_obs)
                # plt.plot(R_icsd[:,0])

    #%% 
    
    
    k_MALM = []
    for i, f in enumerate(selected_files_MALM):
        print('-'*24)
        print(f)
        print('-'*24)

        background_ERT_time = int(f.split('_')[2])
        for j, n in enumerate(ERT_log['Name'][ERT_log['method'] == 'ERT']):
            if 'PRD_ERT_' + str(background_ERT_time) in n:
                index_ERT_backgrd = j
        # ERT_log[ERT_log['Name'].isin(['PRD_ERT_' + str(background_ERT_time)])]
        # ERT_log[ERT_log['Name']==('PRD_ERT_' + str(background_ERT_time))]
        # print(background_ERT_time)
        
        m = 71
        filterSeq8 = args.filter_seq
        if '8returns' in f:
            filterSeq8 = False
            m = None
            
        outMALM = proc.prepare_icsd(imaging_path + inversionPathMALM,
                                        f,
                                        k_indiv_merged[index_ERT_backgrd],
                                        # reprocessed=1, #bool(args.reprocessed),
                                        reprocessed=bool(args.reprocessed),
                                        nVRTe_cols=6*2+1, nVRTe_rows=4*2+1,
                                        filter_seq_rec=args.filter_seq_rec,
                                        filter_seq=filterSeq8,
                                        m=m,
                                        reduce2d=True,
                                        )
        k_MALM.append(outMALM[0])
        [_, imin, R_obs, nodes, R_icsd] = outMALM


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
        plt.close()
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

    icsd_stk = []
    sol_stk = []
    f_names = []
    rms = []
    j = 0

    for i, f in enumerate(selected_files_MALM):
        if (('8returns' in f)):# and (args.filter_seq==True)):
            continue
        
            # icsd, sol, ax, fig = proc.invert_MALM(imaging_path + inversionPathMALM + f,
            #                             wr=args.wr, typ=args.dim, pareto=pareto, show=True,
            #                             prior = prior,
            #                             # fname_sim='VRTeSim_icsd.txt'
            #                             )
        
        print('*'*24)
        print(f)
        print('*'*24)
        pareto = bool(args.pareto)
        prior = False
        # for prior in [True,False]:
        icsd, sol, ax, fig = proc.invert_MALM(imaging_path + inversionPathMALM + f,
                                    wr=args.wr, typ=args.dim, pareto=pareto, show=True,
                                    prior = prior,
                                    # fname_sim='VRTeSim_icsd.txt'
                                    )
        sol_stk.append(sol)
        icsd_stk.append(icsd)
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

        # RMS analysis of the ICSD
        # -------------------------------------------------------------------
        date_title = ERT_log[ERT_log['Name']==f]['datetime'].dt.date.values[0]

        rms.append(
                    icsd.RMSAnalysis(
                                    path=imaging_path + inversionPathMALM + f + '/', 
                                    prefix_name=f,
                                    title = str(date_title)
                                    )
                   )
        
        
        

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
    
    
    # np.savetxt(imaging_path + 
    #             inversionPathMALM + 
    #             'rms_ICSD.txt',
    #             rms)
    
    
    
    fig, ax = plt.subplots(1)
    
    rms_df = pd.DataFrame(rms)
    rms_df.columns = ['rms_ICSD']
    rms_df['datetime'] = ERT_log[ERT_log['Name'].isin(f_names)]['datetime'].values
    rms_df.set_index('datetime',inplace=True)
    rms_df.plot.hist(y='rms_ICSD')
    

    rms_df['rms_ICSD_perc'] = rms_df['rms_ICSD']*100
    rms_df.to_csv(imaging_path + 
                inversionPathMALM + 
                'rms_ICSD.csv'
                )
    
    plt.savefig(imaging_path + 
                inversionPathMALM + 'RMS_ICSD.png', dpi=350)
    
    np.savetxt(imaging_path + 
                inversionPathMALM + 
                'vrte_nodes_in_mesh.txt',
                imin)
    
    #%%
    
    # # if args.filter_seq_rec == False:
    # fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))
    
    # # df_SWC
    
    # ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)

    # surveyPRD.plot_PRD_effect_MALM(k_indiv_merged, df_MALM, ax=ax[1], **kwargs)
    # # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))
    # ax[1].grid('on', which='major', axis='x',color='0.95' )
    # ax[1].grid('on', which='major', axis='y',color='0.95' )
    # plt.xticks(rotation=35)
    # plt.tight_layout()

    # if type(args.cycle) is int:
    #     args.cycle = [args.cycle]
    # plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect' + str([*args.cycle]),
    #             dpi=450)
    # plt.close('all')


    #%%

    proc.plot_scatter_MALM(imaging_path+inversionPathMALM,
                           df_MALM,
                           k_indiv_merged,
                           ERT_log,
                           vmin=-625, vmax=125,
                           f_MALM=f_MALM,
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

    # fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))   
    # ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)
    # surveyPRD.plot_PRD_effect_MALM_diff(k_indiv_merged, df_MALM_diff_b, ax=ax[1])
    # # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))
    # plt.xticks(rotation=35)
    # plt.tight_layout()

    # if type(args.cycle) is int:
    #     args.cycle = [args.cycle]
    # plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect_diff_background' + str([*args.cycle]),
    #             dpi=450)
    # plt.close('all')

    
    
    df_MALM_diff = proc.plot_scatter_MALM_diff(imaging_path+inversionPathMALM,
                                df_MALM,
                                k_indiv_merged,
                                f_MALM,
                                ERT_log,
                                background_diff=True,
                                vmin=25, vmax=150,
                                )
    
    # fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))   
    # ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)
    # surveyPRD.plot_PRD_effect_MALM_diff(k_indiv_merged, df_MALM_diff, ax=ax[1])
    # ax[1].grid('on', which='major', axis='x',color='0.95' )
    # ax[1].grid('on', which='major', axis='y',color='0.95' )
    # # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))
    # plt.xticks(rotation=35)
    # plt.tight_layout()

    # if type(args.cycle) is int:
    #     args.cycle = [args.cycle]
    # plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect_diff' + str([*args.cycle]),
    #             dpi=450)
    # plt.close('all')

    
    
    #%%
    if len(f_MALM) == 2*len(k_indiv_merged):
        df_MALM_diff_StemSoil = proc.plot_scatter_MALM_diff_stem_soil(
                                            imaging_path+inversionPathMALM,
                                            df_MALM,
                                            k_indiv_merged,
                                            f_MALM,
                                            ERT_log,
                                            vmin=25, vmax=150,
                                            )
    plt.close('all')

        

    
    #%%
    dates = ERT_log[ERT_log['Name'].isin(f_names)]['datetime'].dt.strftime("%d %B").values #.strptime.dt.date.values

    def subplot_results_icsd(icsd_stk,f_names,dates,savename,clim=[0,0.1]):
        # if len(f_names)>10:
            # nrows= int(len(f_names)/2)
        # fig , axs = plt.subplots(int(len(f_names)/7),int(len(f_names)/6),#figsize=(18,18),
        #                          sharex=True,sharey=True)
        
        
        nrows = 3
        ncols = 2
        if len(f_names)>6:
            nrows = 3
            ncols = 3
        if len(f_names)>9:
            nrows = 4
            ncols = 4
        if len(f_names)>16:
            ncols = 5
            nrows = 4
        if len(f_names)>20:
            ncols = 5
            nrows = 5
        if len(f_names)>25:
           ncols = 6
           nrows = 5
        if len(f_names)>30:
           ncols = 6
           nrows = 6
           
        fig , axs = plt.subplots(nrows,ncols,
                                 sharex=True,
                                 sharey=True,
                                 #figsize=(5,3)
                                 )

        ax_r= axs.ravel()    
        j = 0          
        for i in range(len(f_names)):
            
            print('oo'*21)
            print(f_names[i])
            print('oo'*21)

            cbarplot = False
            if i == len(f_names)-1:
                cbarplot = True
                
            ax_subi, f_subi = icsd_stk[i].showResults(
                                                        ax=ax_r[j],
                                                        fig_name=str(dates[i]),
                                                        lgd_label='',
                                                        cbarplot=cbarplot,
                                                        clim=[0,0.1],
                                                        vrte_pos=False,
                                                        )
            ax_r[j].set_title(dates[i],fontsize=6)
            ax_r[j].tick_params(axis='y', labelsize=6)

            ax_r[j].set_ylabel('z [m]', fontsize=7)
            if i>0:
                ax_r[j].set_ylabel('')
            if i == len(f_names):
                ax_r[j].set_xlabel('x [m]', fontsize=7)
            else:
                xax = ax_r[j].axes.get_xaxis()
                xax = xax.set_visible(False)
            j += 1 
    
        for i2del in range(len(ax_r)):
            if i2del>=len(icsd_stk):
                fig.delaxes(ax_r[i2del])
                # axes[i, j].remove()
    
        plt.subplots_adjust(wspace=-0.8, hspace=0.8)    
        plt.tight_layout()        
        plt.savefig(imaging_path + 
                    inversionPathMALM + savename, dpi=400,
                    bbox_inches='tight', pad_inches=0
                    )
    

    matching = [i for i, s in enumerate(f_names) if "MALM2" in s]
    matching_bool = [False]*len(f_names)
    matching_bool = np.array(matching_bool)
    matching_bool[matching] = True
    
    subplot_results_icsd(
                         np.array(icsd_stk)[matching_bool],
                         np.array(f_names)[matching_bool],
                         np.array(dates)[matching_bool],
                         clim=[0,0.1],
                         savename='subplots_ICSD.png'
                         )
    subplot_results_icsd(
                        np.array(icsd_stk)[~matching_bool],
                        np.array(f_names)[~matching_bool],
                        np.array(dates)[~matching_bool],
                        clim=[0,0.1],
                        savename='subplots_ICSD_stem.png'
                         )
        
    #%%

    
    # fig, ax = plt.subplots((2), figsize=(8, 5), sharex=(True))   
    # ax0 = surveyPRD.plot_timeline(ERT_log, irr_log, ax=ax[0], dropMALM=True)
    
    # df_MALM_icsd = surveyPRD.plot_PRD_effect_icsd(k_indiv_merged, 
    #                                               imin,
    #                                               df_MALM_icsd, irr_log,
    #                                               ax=ax[1],
    #                                               hours_interval=10,
    #                                               **kwargs)
    # ax[1].grid('on', which='major', axis='x',color='0.95' )
    # ax[1].grid('on', which='major', axis='y',color='0.95' )
    # # ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=5))
    # # ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))
    # plt.xticks(rotation=35)
    # plt.tight_layout()
    
    # # Add a zoom to irrigation
    # # -----------------------------------------------------------------

    # if type(args.cycle) is int:
    #     args.cycle = [args.cycle]
    # plt.savefig(imaging_path + inversionPathMALM + 'PRDeffect_icsd' + str([*args.cycle]),
    #             dpi=450)
    
    
# %%
# import pyvista as pv
# pv.__version__ # '0.34.1

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
    
    detrend = 0.1/4

    inversionPathMALM = surveyPRD.definePaths(args)

    # Read and parse log files
    # -----------------------------------
    f_ERT, f_MALM,  ERT_log, irr_log = surveyPRD.load_log_files(args)

    # test
    # ERT inversion
    # -----------------------------------
    # k_indiv_ERT, df_ERT = run_ERT_ind(args, f_ERT, ERT_log,
    #                                   detrend=detrend)
    

    # TL ERT inversion
    # -----------------------------------
    if args.TL:
        run_TL(f_ERT)

    # %% MALM icsd inversion
    # -----------------------------------
    if args.icsd:
        run_icsd(f_MALM, f_ERT, k_indiv_ERT, inversionPathMALM,
                 detrend=detrend)
