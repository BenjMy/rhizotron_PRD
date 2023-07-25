from scipy.spatial.distance import cdist
from pyPRD import processing as proc
import numpy as np
import pyvista as pv
import argparse
import os
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib as mpl
import surveyPRD


imaging_path = '../imaging_ERT_MALM/'

# %%


def get_cmd():
    parse = argparse.ArgumentParser()

    period = parse.add_argument_group('period')
    period.add_argument('-cycle', '--cycle', nargs='+',
                        help='list of cycle', type=int, default=[-99], required=False)
    period.add_argument(
        '-startD', type=str, help='start date to analyse', default='None', required=False)
    period.add_argument('-endD', type=str, help='end date',
                        default=None, required=False)


    process_param = parse.add_argument_group('process_param')
    process_param.add_argument(
        '-scenario', type=str, help='Scenario', default='4Sbis2', required=False)
    process_param.add_argument(
        '-recErr', type=int, help='Rec. error', default=5, required=False)
    process_param.add_argument(
        '-filter_seq', type=int, help='filter sequence', default=0, required=False)
    process_param.add_argument('-filter_seq_rec', type=int,
                               help='filter sequence rec', default=1, required=False)
    process_param.add_argument(
        '-reprocessed', type=int, help='reprocessed', default=0, required=False)
    process_param.add_argument('-anisotropy', type=float,
                               help='anisotropy procesing', default=1, required=False)
    process_param.add_argument('-pareto', type=int,
                               help='pareto', default=1, required=False)
    process_param.add_argument('-wr', type=float,
                               help='reg. weight', default=1, required=False)
    process_param.add_argument('-obs_err', type=str,
                               help='obs_err. weight', default='const', required=False)
    
    
    args = parse.parse_args()
    return(args)


args = get_cmd()

# cycle 7
# args.startD = '29/6/2022,14:14'
# args.endD = '29/6/2022,15:03'

args.startD = '11/7/2022,15:49'
args.endD = '11/7/2022,17:11' #17h10 8 returns
# args.startD = '8/6/2022,9:59'
# args.endD = '8/6/2022,11:31'


sc_nb = args.scenario



def define_sources_prop(sc_nb):
    
    ls_sc = {
                '1S': {'idS': [59], 'wS': [1], 'rho0':100,'prior': False,'reduce2d':True,
                      },
                '1S_ERT': {'idS': [59], 'wS': [1], 'rho0':None,'prior': False,'reduce2d':True,
                      },
                '2S': {'idS': [59,89], 'wS': [0.5,0.5], 'rho0':100,'reduce2d':True,
                      'prior': False,
                      },
                '3S': {'idS': [59,89,100], 'wS': [1/3,1/3,1/3], 'rho0':100,'reduce2d':True,
                      'prior': False,
                      },
                '4S': {'idS': [59,89,100,110], 'wS': [1/4,1/4,1/4,1/4], 'rho0':100,'reduce2d':True,
                      'prior': False,
                      },
                '4Sbis': {'idS': [59,89,100,110], 'wS': [1/2,(1/2)/3,(1/2)/3,(1/2)/3], 'rho0':100,'reduce2d':True,
                      'prior': False,
                      },
                '4Sbis_ERT': {'idS': [59,89,100,110], 'wS': [1/2,(1/2)/3,(1/2)/3,(1/2)/3], 'rho0':None,'reduce2d':True,
                      'prior': False,
                      },
               '4Sbis2': {'idS': [59,89,100,110], 'wS': [1/10,(9/10)/3,(9/10)/3,(9/10)/3], 'rho0':100,'reduce2d':True,
                     'prior': False,
                     },
        }

    return ls_sc[sc_nb]

scenario = define_sources_prop(sc_nb)

def definePaths(args,sc_nb):
        
    inversionPathMALM = 'inversionMALM_synth/' + sc_nb + '/'
    if args.filter_seq:
        inversionPathMALM = 'inversionMALM_filter_synth/' + sc_nb + '/'
    if args.filter_seq_rec:
        inversionPathMALM = 'inversionMALM_filter_rec_synth/' + sc_nb + '/'
        
    return inversionPathMALM

inversionPathMALM = definePaths(args,sc_nb)

# %%

ERT_log = surveyPRD.load_ERT_survey_log(startDate=args.startD, endDate=args.endD,
                                        cycles=args.cycle)
irr_log = surveyPRD.load_irr_log(startDate=args.startD, endDate=args.endD,
                                       cycles=args.cycle)

selected_files_ERT = ERT_log[ERT_log['method'] == 'ERT']['Name'].to_list()
selected_files_MALM = ERT_log[ERT_log['method'] == 'MALM']['Name'].to_list()



# %%
ax = surveyPRD.plot_timeline(ERT_log, irr_log, dropMALM=True)
# tdelta = (ERT_log['datetime'].iloc[0]-ERT_log['datetime'].iloc[-1])/10
# ax.get_xaxis().set_major_locator(mdates.HourLocator(interval=55))
# ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%d/%m - %Hh"))


# %%
# k.createBatchSurvey(imaging_path + 'filenames_ERT')
# k.invert(parallel=True)

k_indiv_merged = surveyPRD.process_ERT(imaging_path,
                                       selected_files_ERT,
                                       ERT_log,
                                       reprocessed=False, #bool(args.reprocessed),
                                       recip=args.recErr,
                                       )


# df_rms = proc.getR2out(imaging_path,
#                        selected_files_ERT
#                        )


# k_indiv_merged[0].dirname
# k_indiv_merged[0].getR2out()

# df = pd.DataFrame(columns=['name', 'dataset', 'iteration', 'resRMS',
#                            'phaseRMS', 'read', 'rejected', 'success'])

# %%
proc.plot_ERT(k_indiv_merged[0], vmin=-10, vmax=1,
              attr="Sensitivity_map(log10)", index=0)
proc.plot_ERT(k_indiv_merged[0], vmin=0, vmax=50,
              attr="Resistivity(ohm.m)", index=0)

# %%
# k_indiv_merged[0].elec
nodes = k_indiv_merged[0].mesh.node
k_MALM_stk = []
idS = scenario['idS']
wS = scenario['wS']
R_obs_stck = []
f_MALM = []

for f in selected_files_MALM[:]:
    
    print(f)

    a, b, m = [72, 65, 71]
    # test for anisotropy type MALM
    # -----------------------------
    filterSeq = args.filter_seq
    if '8returns' in f:
        filterSeq = False
        m = None
        
    outMALM = proc.prepare_MALM_synth(imaging_path+inversionPathMALM,
                                      f,
                                      k_indiv_merged[0],
                                      reduce2d=scenario['reduce2d'],
                                      nVRTe_cols=6*2+1, nVRTe_rows=4*2+1,
                                      idS=idS,
                                      wS=wS,
                                      filter_seq_rec=args.filter_seq_rec,
                                      filter_seq=filterSeq,
                                      a=a, b=b, m=m,
                                      percent_rec=args.recErr,
                                      rho0= scenario['rho0'],
                                      )
    if type(outMALM) == tuple:
        [k_MALM, imin, R_obs, nodes, R_icsd, R_sim, idS, wS] = outMALM
    else:
        [k_MALM, imin, _, nodes, _, _, _, _]  = outMALM
    mesh = k_MALM.mesh.df
    
    
    # plt.plot(R_obs)
    # plt.plot(R_sim[:,0])
    
    # np.shape(R_sim)
    
    
    # k_MALM.param['node_elec']
    #%%
    imin, nodes, grid = proc.create_source_grid(mesh, k_MALM,
                                                nVRTe_cols=6*2+1, nVRTe_rows=4*2+1)
    
    # fig, ax = plt.subplots(1,1)
    # ax.scatter(grid[:,0],grid[:,2])
    # ax.set_xlim([0,0.47])
    # ax.set_ylim([0,0.5])
    # ax.set_xlabel('x (m)')
    # ax.set_xlabel('z (m)')
    # elecs = k_indiv_merged[0].surveys[0].elec
    # ax.scatter(elecs['x'],elecs['z'])
    # ax.set_xlim([0,0.47])
    # ax.set_ylim([0,0.50])
    # plt.show()
    
    
    text_file = open("vrte_positions.txt", "wt")
    for i, g in enumerate(grid):
        pt_nb = i + 81
        defaut_line = 'Point( {} ) = {{{}, {}, {}, cl1}}; \n'.format(pt_nb, g[0],g[1],g[2])
        n = text_file.write(defaut_line)
    text_file.close()

    text_file2 = open("vrte_in_mesh_str.txt", "wt")
    stri = ''
    for i in range(len(grid)+81):
        # if (i>90+72 & i<90+81):
        print(i)
        if i==(90+72):
            stri += str(i+1)
        else:
            stri += str(i+1) + ','
    n = text_file2.write(stri)
    text_file2.close()



    k_MALM_stk.append(k_MALM)
    xyz_sol = nodes[imin[idS]]
    
    
    # xyz_sol = nodes[imin[7*9-4]]
    # xyz_sol = nodes[imin[7*9+26]]
    # xyz_sol = nodes[imin[100]]
    # xyz_sol = nodes[imin[110]]

    # xyz_sol = nodes[imin[80]]
    # xyz_sol = nodes[imin[82]]
    # xyz_sol = nodes[imin[84]]


    # xyz_sol = nodes[imin[idS]]
    
    # if len(R_obs_stck)>0:
    #     if len(R_obs)==len(R_obs_stck[0]):
    #         R_obs_stck.append(R_obs)
    #         f_MALM.append(f)
    # else:
    #     R_obs_stck.append(R_obs)
    #     f_MALM.append(f)


    # test for anisotropy type MALM
    # -----------------------------
    # if '8returns' in f:
    #     df_anisotropy = surveyPRD.anisotropy(k_MALM)
    #     # plot_ABMN_pos(k_MALM,df_anisotropy)
    #     df_anisotropy_new = surveyPRD.normalised_ratio(k_MALM,df_anisotropy)
    #     eucl_dist_A_sol, angle_A_sol = surveyPRD.sol_polar(xyz_sol)
    #     surveyPRD.plot_anisotropy_polar(df_anisotropy, eucl_dist_A_sol, angle_A_sol,
    #                                     imaging_path+inversionPathMALM, norm=True)
    #     surveyPRD.plot_anisotropy_polar(df_anisotropy, eucl_dist_A_sol, angle_A_sol,
    #                                     imaging_path+inversionPathMALM,  norm=False)




# for ff in f_MALM:
#     if'return' in ff:
#         f_MALM.remove(ff)
            
# df_MALM = surveyPRD.add2df_MALM(R_obs_stck,
#                                 f_MALM,
#                                 ERT_log)

# proc.plot_scatter_MALM(imaging_path+inversionPathMALM,
#                         df_MALM,
#                         k_indiv_merged,
#                         ERT_log,
#                         vmin=-125, vmax=550,
#                         )

# k_MALM_stk[0].sequence = k_MALM_stk[0].surveys[0].df[['a','b','m','n']]


# proc.plot_scatter_MALM(imaging_path+inversionPathMALM,
#                        df_MALM,
#                        k_indiv_merged,
#                        ERT_log,
#                        vmin=-125, vmax=550,
#                        )


# k_MALM.df


        
        
# %%
dim = '3d'
if scenario['reduce2d']:
    dim = '2d'
         
for i, f in enumerate(selected_files_MALM[:]):
    f = selected_files_MALM[0]
    # i=0
    m0_MALM = proc.m0_MALM(imaging_path + inversionPathMALM + f,
                            method_m0='F1', typ=dim
                            )
    proc.plot_m0_MALM(k_MALM_stk[i], m0_MALM, show=False)
    # print(xyz_sol)

# %%
    # m0_MALM = proc.m0_MALM(imaging_path + inversionPathMALM + f,
    #                         method_m0='Pearson', typ=dim
    #                         )
    # proc.plot_m0_MALM(k_MALM_stk[i],m0_MALM,show=False)
    # print(xyz_sol)

# %%
for i, f in enumerate(selected_files_MALM[:]):
    print(xyz_sol)

    # f = selected_files_MALM[2]
    pareto = False
    sol = proc.invert_MALM(imaging_path + inversionPathMALM + f,
                           pareto=args.pareto, wr=args.wr, typ=dim,
                           prior=scenario['prior'],
                           obs_err=args.obs_err
                           )
    
    if args.pareto:
        proc.plot_MALM(k_MALM_stk[i], sol[1], show=False, ext=['png'])
        # proc.plot_MALM(k_MALM[i],sol,show=True,ext=['svg'])
    else:
        proc.plot_MALM(k_MALM_stk[i], sol[1], show=False, ext=['png'],
                       )
        # proc.plot_MALM(k_MALM_stk[i], sol, show=False, ext=['png'],
        #                log_scale=True)
        # proc.plot_MALM(k_MALM[i],sol,show=True,ext=['svg'])
      
    # f, axs = plt.subplots(1)
    ax, fig = sol[0].showResults(
                                # ax=axs,
                                # fig_name=str(dates[i]),
                                lgd_label='',
                                # cbarplot=cbarplot,
                                clim=[0,1],
                                vrte_pos=False,
                                )
    ax.scatter(xyz_sol[:,0],xyz_sol[:,2],color='c',marker='+',s=100)
    fig.savefig(imaging_path + inversionPathMALM + f + 'ICSD.png')
    

if bool(args.anisotropy):
    k_MALM_anisotropy = []
    m0_anisotropy = []
    for i, f in enumerate(selected_files_MALM):
        if '8returns' in f:
            # b_return = [1,4,8,17,20,24,33,36,40,49,52,64]
            b_return = [1,8,17,24,36,49,64]
            for k in range(len(b_return)):
                print('*'*24)
                print(k,f)
                print('*'*24)
                fnew = f.split('.csv')[0] + '_B'+ str(b_return[k])
                background_ERT_time = int(f.split('_')[2])
                for j, n in enumerate(ERT_log['Name'][ERT_log['method'] == 'ERT']):
                    if 'PRD_ERT_' + str(background_ERT_time) in n:
                        index_ERT_backgrd = j
                outMALM = proc.prepare_MALM_synth(imaging_path + inversionPathMALM,
                                                    fnew,
                                                    k_indiv_merged[index_ERT_backgrd],
                                                    nVRTe_cols=6*2+1, nVRTe_rows=4*2+1,
                                                    idS=idS,
                                                    wS=wS,
                                                    filter_seq_rec=args.filter_seq_rec,
                                                    filter_seq=args.filter_seq,
                                                    reduce2d=True,
                                                    percent_rec=args.recErr,
                                                    a=72, b=b_return[k], m=None,
                                                    rho0= scenario['rho0'],
                                                    retElec=k_indiv_merged[0].elec[['x','z']].iloc[b_return[k]-1].to_numpy(),
                                                    )
                k_MALM_anisotropy.append(outMALM[0])
                
                m0i,ax,fig = proc.m0_MALM(imaging_path + inversionPathMALM + fnew,
                                       method_m0='F1', typ='2d',
                                       show=True
                                       )
                ax.scatter(xyz_sol[:,0],xyz_sol[:,2],color='c',marker='+',s=100)
                # plt.close()
                fig.savefig(os.path.join(imaging_path,inversionPathMALM, fnew,'m0' + '_breturn' + str(b_return[k]) + '.png')
                            , dpi=300)
                m0_anisotropy.append(m0i)
                
                icsd, sol, ax, fig = proc.invert_MALM(imaging_path + inversionPathMALM + fnew,
                                            wr=args.wr, typ='2d', pareto=args.pareto, show=True,
                                            prior = False,                                                        
                                            clim=[0,1],
                                            retElec=k_indiv_merged[0].elec[['x','z']].iloc[b_return[k]-1].to_numpy()
                                            # fname_sim='VRTeSim_icsd.txt'
                                            )
                
                ax.scatter(xyz_sol[:,0],xyz_sol[:,2],color='c',marker='+',s=100)
                fig.savefig(os.path.join(imaging_path,
                                         inversionPathMALM, 
                                         fnew,
                                         '_breturn' + str(b_return[k]) + 
                                         '_icsd.png'), 
                            dpi=300)
                plt.close('all')
            plt.close('all')
        
        
# if pareto:
    # proc.plot_MALM(k_MALM_stk[i],sol[0])
    # else:
    # proc.plot_MALM(k_MALM_stk[i], sol,
    #                 )
    # print(xyz_sol)

