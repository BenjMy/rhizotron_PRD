# %%capture
# !pip install numpy GitPython pyvistaqt itkwidgets 
import os
import numpy as np

# import pygimli as pg

from pyCATHY import cathy_tools 
# 
from pyCATHY.plotters import cathy_plots as cplt 
from pyCATHY.importers import cathy_inputs as  in_CT
from pyCATHY import cathy_utils as utils

from pyCATHY import meshtools as mt 

# from pyCATHY import cathy_utils
# from pyCATHY import rhizo_tools
from pyCATHY import rhizo_tools
rhizo = rhizo_tools.rhizotron()
# 
from pyCATHY.DA import cathy_DA
cathyDA = cathy_DA.DA()

import pyvista as pv

# https://gitlab.kitware.com/paraview/paraview/-/issues/20019


# RMSE is average differences over the observations
# RMSE_ERT (mean of all observations=all ert mesh nodes), RMSE_SWC, need to differentiate between observations
import argparse
from rhizo_scenarii.DA_scenarii import load_scenario 


def get_cmd():
    parser = argparse.ArgumentParser(description='Process')
    # parser.add_argument('sc', type=int, help='scenario nb')
    #  #Feddes_irr, ref_scenario_irr_het_soil_f5, hetsoil_irr, Archie, Feddes_RWU
    # -------------------------------------------------------------------------------------------------------
    parser.add_argument('-study','--study', type=str, help='study selection', required=False, default='test')
    parser.add_argument('-sc','--sc', help='scenario nb', required=False, default=3)
    parser.add_argument('-nens','--nens', type=int, help='nb of ensemble', required=False, default=32)
    parser.add_argument('-openLoop','--openLoop', type=int, help='openLoop',default=0)
    parser.add_argument('-freq','--freq', type=int, help='freq',default=1)
    parser.add_argument('-DAtype',type=str, help='type of DA',default='enkf_analysis_inflation')
    parser.add_argument('-damping',type=float, help='damping factor',default=1)
    parser.add_argument('-dataErr',type=float, help='error data',default=0.01)
    parser.add_argument('-parallel',type=int, help='parallel computing',default=1)
    # parser.add_argument('-DAtype',type=str, help='Type of DA',default='enkf_analysis_inflation')
    # parser.add_argument('-DAtype',type=str, help='Type of DA',default='pf')
    
    
    

    args = parser.parse_args()

    return(args)


if __name__ == '__main__':
    
    args = get_cmd()


    #%% Scenario description
    
    """We initiate a CATHY object; if the the CATHY src files are not included within the 'path2prj', they are automatically fetched from the gitbucket (providing the notebook is initiate locally with an internet connection)"""
    
    path2prj ='rhizo_models'    
    
    scenarii = load_scenario(study=args.study)
    
    Snb = int(args.sc) # 3 # -1
    NENS = int(args.nens) #128 128*2
    open_loop_run = bool(args.openLoop)
    parallel = bool(args.parallel)
    DA_type = args.DAtype
    prj_name = list(scenarii)[Snb] 
    
    

    
    
    
    #%% load solution parameters
    
    import pickle
    with open(os.path.join(path2prj, 
                           scenarii[prj_name]['ref_scenario'], 
                           scenarii[prj_name]['ref_scenario'] + '.pkl'), 'rb') as f:
        solution_parm = pickle.load(f)
    print(solution_parm)
    
    #%% Update input files
    pmin = solution_parm['pmin']
    
    # os.getcwd()
    nb_of_days = solution_parm['nb_of_days']
    """We initiate a CATHY object; if the the CATHY src files are not included within the 'path2prj', they are automatically fetched from the gitbucket (If the notebook is initiate locally with an internet connection)"""
    
    path2prj ='mary_rhizo_withDA'
    times = list(np.arange(0,solution_parm['nb_of_days']*86400,solution_parm['delta_t']))
    simu_time_max = times[-1]+1
    
    #%% ERT settings
    
    pathERT =  solution_parm['pathERT_results'] # '/home/ben/Documents/CATHY/pyCATHY/ERT_fwd_DA_sol_rhizo/'
    path_inputs_rhizo = '/home/ben/Documents/CATHY/pyCATHY/inputs_rhizo/'
    prjERT = scenarii[prj_name]['ref_scenario']
    
    meshERT = path_inputs_rhizo + 'mesh/BaseRhizo.vtk'
    sequenceERT = path_inputs_rhizo +  'sequences/SequenceERT_Rhizo_71Elecs.shm'
    data_format = 'pygimli'
    elecs= np.genfromtxt(path_inputs_rhizo + 'mesh/elecsXYZ.csv', delimiter=",",skip_header=1)
    
    
    resample = np.arange(0,len(times),args.freq)   
    time_sampled = [times[i] for i in resample]


    
    #%% Create new project dir to save DA results
    
    
    prj_name = (list(scenarii)[Snb] +
              '_F' + str(args.freq) +
              '_OL' + str(args.openLoop) +
              '_NENS' + str(args.nens) +
              '_DAtyp' + args.DAtype + 
              '_alpha' + str(args.damping) +
              '_DataErr' + str(args.dataErr) 
              )
    
    if 'inflation' in args.DAtype:
        prj_name += '_alpha' + str(args.damping)
    
   
    path2prj ='rhizo_withDA'
    simu_DA = cathy_tools.CATHY(dirName=path2prj, 
                                prj_name=prj_name,
                                notebook=False) #
    
    
    #%% read scale data observation
    
    rhizo.prepare_atmbc_from_weight(simu_DA)
    
    weight    = {
                'data_type': 'weight', # units
                'units': '$\kg$', # 
                'instrument': '', # 
                'data_format': ''
                  }


    # for i, tt in enumerate(time_sampled):
        

    #     filename = os.path.join(pathERT, prjERT, 'ER_predicted_sol_t' + str(i) + '.csv')
            
    #     # need to call read_real_data as many times as variable to perturbate
    #     # return a dict merging all variable perturbate to parse into prepare_DA
    #     data_measure = simu_DA.read_observations(filename=filename, 
    #                                             data_type = 'ERT', 
    #                                             data_err = args.dataErr, # instrumental error
    #                                             show=True,
    #                                             tA=tt,
    #                                             obs_cov_type='data_err', #data_err
    #                                             elecs=elecs,
    #                                             meta=ERT
    #                                             ) # data_err  reciprocal_err

    #%% read ERT data observation
    
    ERT    = {
                'data_type': '$ERT$', # units
                'units': '$\Ohm$', # units transfer_resistances
                'forward_mesh_vtk_file': meshERT, # units transfer_resistances
                'sequenceERT': sequenceERT, # units transfer_resistances
                'instrument': 'Syscal', # units transfer_resistances
                'data_format': data_format, # units transfer_resistances
                  }
        
    for i, tt in enumerate(time_sampled):
        
        try:
            # csv file observation generated by resipy
            filename = os.path.join(pathERT, prjERT, 'ER_predicted_sol_t' + str(i) + '.csv')
                
            # need to call read_real_data as many times as variable to perturbate
            # return a dict merging all variable perturbate to parse into prepare_DA
            data_measure = simu_DA.read_observations(filename=filename, 
                                                    data_type = 'ERT', 
                                                    data_err = args.dataErr, # instrumental error
                                                    show=True,
                                                    tA=tt,
                                                    obs_cov_type='data_err', #data_err
                                                    elecs=elecs,
                                                    meta=ERT) # data_err  reciprocal_err
        except:
            pass
        # use reciprocal errors to compute the covariance matrice
        # on diagonal --> 1/abs(dict_obs['data']['recipError']
        # len(data_measure)

   
    #%% Update parameters


    # simu_DA = update_rhizo_inputs(simu_DA, 
    #                               nb_of_days=nb_of_days,
    #                               solution = solution_parm,
    #                               tobs=time_sampled) # ZROOT = 0.2
    rhizo.set_defaults(simu_DA)
    rhizo.update_atmbc_PRD(simu_DA)

    #%% perturbated variable 
    

    parm_per = rhizo_tools.perturbate_rhizo(cathyDA,simu_DA,
                                            scenarii[list(scenarii)[Snb]],
                                            prj_name=prj_name,
                                            NENS=args.nens)
    
    # sampling_type
    # import numpy as np
    # parm_sampling = np.random.normal(0.001,scale=0.25, size=32)
    # parm_sampling = np.random.normal(0.25,scale=0.25, size=32)

    #%%
    
    if not hasattr(solution_parm,'rFluid_Archie'):
        simu_DA._mapping_petro_init()
        solution_parm.update(simu_DA.Archie_parms)
           
    simu_DA.set_Archie_parm(
                            rFluid_Archie=solution_parm['rFluid_Archie'],
                            a_Archie=solution_parm['a_Archie'],
                            m_Archie=solution_parm['m_Archie'],
                            n_Archie=solution_parm['n_Archie']
                            )
    
    #%% Run assimilation cycle
    
    simu_DA.run_processor(DAFLAG=1,
                          DTMIN=10,DELTAT=10,
                          parallel=parallel,    
                          dict_obs= data_measure,
                          list_assimilated_obs = 'all', # default
                          list_update_parm= scenarii[list(scenarii)[Snb]]['listUpdateParm'],
                          DA_type = DA_type, #'pf_analysis', # default
                          dict_parm_pert=parm_per,
                          open_loop_run=open_loop_run,
                          threshold_rejected=80,
                          damping=args.damping                    
                          )