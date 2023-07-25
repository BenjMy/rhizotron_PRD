import warnings
warnings.filterwarnings('ignore')
import os
import sys
sys.path.append((os.path.relpath('../src'))) # add here the relative path of the API folder
testdir = './rawData/'
from resipy import Project
import numpy as np
import pyvista as pv
pv.set_plot_theme("document")
# pv.global_theme.colorbar_horizontal.width = 0.2  
pv.global_theme.colorbar_orientation = 'vertical'  
pv.global_theme.font.label_size = 40  
pv.global_theme.font.size = 48

# pv.global_theme.colorbar_vertical.width = 0.45 
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib as mpl

import multiprocessing
from multiprocessing import Process, Manager
import functools
from functools import partial
import pandas as pd
import math
import pickle
import pathlib
from icsd3d.icsd3d_class import iCSD3d as i3d 

MESH = "mesh/mesh_rhizo_resipy_vrte_v7.msh"


def R2(y_fit,y,x=False):
    # residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)
    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    
    # if x:
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(y_fit,y)
    return r_value, p_value

    # else:
    #     return r2



def prepare_MALM(path, filename, k_ERT, parallel=False,
                 reduce2d=False,
                 **kwargs):
    
    
    if not os.path.exists(path + filename):
        os.makedirs(path + filename)

    k_MALM, R_obs = create_MALM_project(path,filename,**kwargs)

    return R_obs


def _prepare_MALM_(path, filename, k_ERT, parallel=False,
                 reduce2d=False,reprocessed=False,
                 **kwargs):
    
    if not os.path.exists(path + filename):
        os.makedirs(path + filename)

    k_MALM, R_obs = create_MALM_project(path,filename,**kwargs)

    mesh = k_MALM.mesh.df
    imin, nodes, grid = create_source_grid(mesh, k_MALM, **kwargs)
    k_MALM = set_MALM_root_seq(k_MALM)
    
    
    if parallel:
        R_icsd = par_fwd_vrte(k_MALM, k_ERT, imin)
    else:
        R_icsd, R_sim = fwd_vrte(k_MALM, k_ERT, imin)
    
    # np.savetxt(path + 'inversionMALM/' + filename + '/VRTeSim_sim.txt', np.hstack(R_sim), fmt='%1.2e')
    # np.savetxt(path + 'inversionMALM/' + filename + '/VRTeSim_icsd.txt', np.hstack(R_icsd), fmt='%1.2e')

    if reduce2d:
        export4icsd(path + filename,
                R_icsd,R_obs,imin, nodes,reduce2d)
    else:
        export4icsd(path + filename,
                R_icsd,R_obs,imin, nodes,reduce2d)
        
    return k_MALM, imin, R_obs, nodes, R_icsd

def prepare_icsd(path, filename, k_ERT, parallel=False,
                 reduce2d=False,reprocessed=False,
                 **kwargs):
    ''' Prepare icsd files '''
    
    

    if os.path.exists(os.path.join(path, 
                                   filename, 
                                   'VRTeCoord.txt')):
        if reprocessed:
            outMALM = _prepare_MALM_(
                                        path,
                                        filename,
                                        k_ERT,
                                        reduce2d=reduce2d,
                                        **kwargs,
                                    )
            [k_MALM, imin, R_obs, nodes, R_icsd] = outMALM
            k_MALM.sequence = k_MALM.sequence.to_numpy()
            
            # if 'b' in kwargs:
                # k_MALM.saveProject(path + filename + '_return' + str(kwargs['b']) + 'backup')
            # else:
            k_MALM.saveProject(path + filename + 'backup')   

        else:
            print('reload file from *.resipy backup')
            k_MALM = Project(dirname=os.path.join(path, filename), typ='R3t')
            
            # if 'b' in kwargs:
            #     k_MALM.loadProject(os.path.join(path, filename  + '_return' + str(kwargs['b']) + "backup.resipy"))
            # else:
            k_MALM.loadProject(os.path.join(path, filename  + "backup.resipy"))

            imin, nodes, grid = create_source_grid(k_MALM.mesh.df, k_MALM, **kwargs)
            
            outMALM = [k_MALM, imin, None , nodes, None]

    else:
            outMALM = _prepare_MALM_(
                                        path,
                                        filename,
                                        k_ERT,
                                        reduce2d=reduce2d,
                                        **kwargs,
                                    )
            [k_MALM, imin, R_obs, nodes, R_icsd] = outMALM
            k_MALM.sequence = k_MALM.sequence.to_numpy()
            k_MALM.saveProject(path + filename + 'backup')   

        
    return outMALM



def create_ERT_survey_pg(pathERT,sequence,mesh,**kwargs):

    
    import pygimli as pg
    from pygimli.physics import ert
    import pygimli.meshtools as mt


    noise_level = 5
    if 'noise_level' in kwargs:
        noise_level = kwargs['noise_level']
    
    print('**'*10)
    print(noise_level)
    
    
    # isExist = os.path.exists(pathERT)
    # if not isExist:
    #   # Create a new directory because it does not exist 
    #   os.makedirs(pathERT) 
    # fname_mesh = './ERT_fwd_DA_sol_rhizo/test_pg/BaseRhizo_Vrte.msh'
    # mesh3d=  mt.readGmsh(mesh, verbose=False)
    # pre, ext = os.path.splitext(mesh)
    # print(os.rename(mesh, pre + '.msh'))
    pre, ext = os.path.splitext(mesh)
    try:
        mesh3d =  mt.readGmsh(pre + '.msh', verbose=True)
    except:
        try:
            mesh3d =  pg.load(pre + '.bms', verbose=True)
        except:
            raise ValueError('Cannot read '+ mesh + ': valid extensions are .msh and .bms')

    # fname_seq = '/ERT_fwd_DA_sol_rhizo/test_pg/SequenceERT_Rhizo_72Elecs.shm'
    # shm = pg.load(fname_seq)
    
    if 'shm' in sequence:
        scheme = pg.load(sequence)
    elif 'txt' in sequence:
        shm = pd.read_csv(sequence, delimiter=' ', header=None)
        shm = shm.to_numpy()
        shm_new = np.subtract(shm,1)
        scheme = pg.DataContainerERT()
        if 'dict_ERT' in kwargs:
            scheme.setSensorPositions(kwargs['dict_ERT']['elecs'])
        for i, elec in enumerate("abmn"):
                scheme[elec] = shm_new[:,i]
    else:
        raise ValueError('Sequence file format not recognized use .shm or .txt as a list of quadrupoles')

    # scheme["k"] = ert.createGeometricFactors(scheme)
    res0 = 1
    if 'res0' in kwargs:
        res0 = kwargs['res0']

    # het = ert.simulate(mesh3d, res=res0, scheme=shm, sr=False, noiseLevel=5,
    #                    calcOnly=True, verbose=True)
    # het.set('k', 1.0/ (hom('u') / hom('i')))
    # het.set('rhoa', het('k') * het('u') / het('i'))
    # het.save('simulated.dat', 'a b m n rhoa k u i')
        
    if len(res0) != len(mesh3d.cells()):
        raise ValueError('wrong initial resistivity input')
    
    het = ert.simulate(mesh3d, res=res0, scheme=scheme, 
                        calcOnly=False, verbose=False, 
                        noiseLevel=noise_level)
    # het = ert.simulate(mesh3d, res=np.ones(len(mesh3d.cells())), scheme=scheme, 
    #                    calcOnly=False, verbose=False, 
    #                    noiseLevel=noise_level)
    
    # het['rhoa']

    return het



def prepare_MALM_synth(path, filename, k_ERT, idS=[10], parallel=False,
                       reduce2d=False, rho0=None,
                       **kwargs):
    
    if os.path.exists(os.path.join(path, 
                                   filename, 
                                   'VRTeCoord.txt')):
        k_MALM = Project(dirname=os.path.join(path, filename), typ='R3t')
        k_MALM.loadProject(os.path.join(path, filename  + "backup.resipy"))
        outMALM = [k_MALM, None, None , None, None, None, None, None]
        
        return outMALM

    else:
        k_MALM, _ = create_MALM_project(path,filename,**kwargs)
        
        
        mesh = k_MALM.mesh.df
        imin, nodes, grid = create_source_grid(mesh, k_MALM, **kwargs)
        
        k_MALM = set_MALM_root_seq(k_MALM)
        # k_MALM.sequence
        
        if rho0:
            k_MALM.mesh.df['res0'] = rho0
        else:
            k_MALM.mesh.df['res0'] = k_ERT.meshResults[0].df['Resistivity(ohm.m)']
    

        # adding multiple_sources
        def multiple_sources(k,idS=[],wS=None):
            seqMALM =  k.surveys[0].df[['a','b','m','n']]
            k.param['node_elec']
            RobsSi = np.zeros((seqMALM.shape[0], len(idS)))
            if len(idS)>0:
                for i, s in enumerate(idS):
                    k.param['node_elec'][1][71] = imin[s] # the source to be found (CHANGE HERE)
                    # k.param['node_elec'][1][int(seqMALM['a'][i])-1] = imin[s] # the source to be found (CHANGE HERE)
                    # k.param['node_elec'][1][71] = s # the source to be found (CHANGE HERE)
                    # k.mesh['res0']
                    k.surveys[0].df[['a','b','m','n']]
                    k.forward()
                    RobsSi[:,i] = k.surveys[0].df['resist'].values
                    # RobsSi[:,i] = 100
                for i, s in enumerate(idS):
                    if wS is None:
                        RobsSi[:,i] = RobsSi[:,i]/len(idS)
                    else:
                        RobsSi[:,i] = RobsSi[:,i]*wS[i]
            Robs = np.sum(RobsSi,axis=1)
            return Robs
        
        # adding noise
        def addnoise(x, level=0.05):
            return x + np.random.randn(1)*x*level
        
        wS=None
        if 'wS' in kwargs:
            wS = kwargs['wS']
        if 'idS' in kwargs:
            idS = kwargs['idS']
        
        k_MALM.param['node_elec']

        # len(k_MALM.sequence)
        R_obs = multiple_sources(k_MALM,idS,wS)
        # addnoise = np.vectorize(addnoise)
        R_obs = addnoise(R_obs, 0.01)    # R_obs = addnoise(R_obs, 0.01)          
        
        if parallel:
            R_icsd = par_fwd_vrte(k_MALM, k_ERT, imin)
        else:
            R_icsd, R_sim = fwd_vrte(k_MALM, k_ERT, imin, rho0)
        


        if reduce2d:
            export4icsd(path + filename,
                    R_icsd,R_obs,imin, nodes,reduce2d)
        else:
            export4icsd(path + filename,
                    R_icsd,R_obs,imin, nodes,reduce2d)
    
        
        k_MALM.sequence = k_MALM.sequence.to_numpy()
        k_MALM.saveProject(path + filename + 'backup')   
    
        # k_MALM.saveProject(path)   
    
            
        return k_MALM, imin, R_obs, nodes, R_icsd, R_sim, idS, wS



def plot_scatter3d_onmesh(point_cloud,interpolated,scalars='csd',pl=None,**kwargs):
    
    # if pl is None:
    #     pl = pv.Plotter()
    
    log_scale = False
    if 'log_scale' in kwargs:
        log_scale = kwargs['log_scale']
        
    pl.add_mesh(point_cloud, scalars=scalars, point_size=20.0, render_points_as_spheres=True,
                log_scale=log_scale)
    pl.add_mesh(interpolated, scalars=scalars,show_edges=True, lighting=False,
                log_scale=log_scale)
    
    actor = pl.show_bounds(grid='front', location='outer',
                            all_edges=True)
    
    # pl.show()
    return pl

def get_vrte_pos(k):
    path_coords = k.dirname + '/..'
    vrte_coords = pd.read_csv(path_coords + '/VRTeCoord.txt',sep='\t', header=None)
    points =  np.c_[vrte_coords[0],-np.ones(len(vrte_coords[0]))*0.005,vrte_coords[1]]
    point_cloud = pv.wrap(points)
    return point_cloud
    
def screenshot_fig_MALM(k_MALM,pl,ext,attr,index, **kwargs):
    pl.view_xz(negative=True)
    #actor = pl.show_grid()
    actor = pl.show_bounds(grid='back', location='front',
                            all_edges=True, color='black',use_2d=True)
    # pl.remove_scalar_bar()
    # pl.add_scalar_bar()
    
    log_scale = ''
    if 'log_scale' in kwargs:
        if kwargs['log_scale']:
            log_scale = '_log'

    prior = ''
    if 'prior' in kwargs:
        if kwargs['prior']:
            log_scale = '_prior'
            
    for e in ext:
        if 'svg' in e:
            pl.save_graphic(k_MALM.dirname + '/' + attr + log_scale + prior + '_t' + str(index) + ".svg")  
            # pl.save_graphic(k_MALM.dirname + '/' + attr + '_t' + str(index) + ".eps")  
            # pl.show(screenshot=k_MALM.dirname + '/' + attr + '_t' + str(index) + ".png")
        if 'png' in e:
            pl.show(screenshot=k_MALM.dirname + '/' + attr + log_scale + prior + '_t' + str(index) + ".png")
            # pl.screenshot(k_MALM.dirname + '/' + attr + '_t' + str(index) + ".png")
    # pl.window_size = window_size
    
def plot_MALM(k_MALM,sol,pl=None,ext='png',show=False, index=0, **kwargs):
    
    pv.set_plot_theme("document")
    off_screen = True
           
    if show:
        off_screen=False
    if pl is None:
        pl = pv.Plotter(off_screen=off_screen,window_size=([2048, 1536]))
        
    pvmesh = pv.read(k_MALM.dirname+'/fwd/forward_model.vtk') # read in temporary file 
    # points =  np.c_[nodes[imin,0],-np.ones(len(nodes[imin,2]))*0.005,nodes[imin,2]]
    point_cloud = get_vrte_pos(k_MALM)
    attr = 'csd'
    
    point_cloud[attr] = sol.x
    # if log_scale:
    #     point_cloud[attr] = np.log(sol.x)
   
    interpolated = pvmesh.interpolate(point_cloud, radius=0.05)
    plot_scatter3d_onmesh(point_cloud,interpolated,scalars=attr,pl=pl,
                          **kwargs)
        
    if show:
        # pl.show()
        screenshot_fig_MALM(k_MALM,pl,ext,attr,index,**kwargs)
    else:
        screenshot_fig_MALM(k_MALM,pl,ext,attr,index,**kwargs)
        return pl
    
def plot_m0_MALM(k_MALM,m0,pl=None,ext=['png'],show=True, index=0,
                 **kwargs):   
    
    pv.set_plot_theme("document")

    off_screen = True
    if show:
        off_screen=False
    if pl is None:
        pl = pv.Plotter(off_screen=off_screen,window_size=([2048, 1536]))
        
    pvmesh = pv.read(k_MALM.dirname+'/fwd/forward_model.vtk') # read in temporary file 
    point_cloud = get_vrte_pos(k_MALM)
    attr = 'm0'
    # len(m0)
    point_cloud[attr] = m0
    
    # if log_scale:
    #     point_cloud[attr] = np.log(m0)
        
    interpolated = pvmesh.interpolate(point_cloud, radius=0.05)

    pl = plot_scatter3d_onmesh(point_cloud,interpolated,scalars=attr,pl=pl)
    pl.show_bounds(all_edges=True)

    if show:
        # pl.show()
        screenshot_fig_MALM(k_MALM,pl,ext,attr,index)
    else:
        screenshot_fig_MALM(k_MALM,pl,ext,attr,index)
        return pl
    

def m0_MALM(path,method_m0='F1', show=False, **kwargs): 
    
    typ = '2d'
    if 'typ' in kwargs:
        typ = kwargs['typ']
    
    icsd=i3d(dirName=path +'/')   
    icsd.regMesh='unstrc'
    icsd.type=typ
    icsd.createSurvey()
    
    
    if show:
        m0, ax, f = icsd.estimateM0(method_m0=method_m0,show=show)
        return m0, ax, f
    else:
        m0 = icsd.estimateM0(method_m0=method_m0,show=show)
        return m0
        
    
    
    # icsd=i3d(dirName=path2files)   
    # icsd.type='2d'
    # icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
    # icsd.wr=1 #weight regularization
    # icsd.alphaSxy=True
    # icsd.x0_prior=True
    # icsd.x0_ini_guess=False # initial guess
    # #icsd.icsd_init()

    # icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim_sim.txt')
    # # icsd.plotElecs=False

    # m0 = icsd.estimateM0(method_m0='F1',show=True)
    
    
    # return m0

def invert_MALM(path,pareto=False,show=False,**kwargs): 
    
    typ = '3d'
    wr = 1
    obs_err = 'sqrt'
    if 'wr' in kwargs:
        wr = kwargs['wr']
    if 'typ' in kwargs:
        typ = kwargs['typ']
    if 'prior' in kwargs:
        prior = kwargs['prior']
    if 'obs_err' in kwargs:
        obs_err = kwargs['obs_err']
        
    icsd=i3d(dirName=path +'/')   
    icsd.regMesh='unstrc'
    icsd.type=typ
    # #icsd.mesh='invdir/fwd/forward_model.vtk'
    icsd.obs_err=obs_err # sqrt choose between constant weight and w = 1/sqrt(abs(obs))
    icsd.wr=wr #weight regularization
    # icsd.alphaSxy=False
    icsd.x0_prior=prior
    icsd.x0_ini_guess=prior # initial guess
    icsd.createSurvey()
    icsd.pareto_MinErr=1e-2
    icsd.pareto_MaxErr=1e2
    # icsd.pareto_weights = np.array([x*(10**y) for y in range(-2,3)  for x in range(1,10)])[:-8] 

    # sol = icsd.invert(pareto=pareto,**kwargs) 
    
    if show:
        sol, ax, f = icsd.invert(pareto=pareto,show=show,**kwargs) 
        return icsd, sol, ax, f
    else:
        sol = icsd.invert(pareto=pareto,**kwargs) 
        return icsd, sol
        
    
    
    # return sol

def export4icsd(path,R_sim,R_obs,imin,nodes,reduce2d=False):
    # export for iscd dev analysis
    np.savetxt(path + '/VRTeSim.txt', np.hstack(R_sim), fmt='%1.5e')
    np.savetxt(path  +'/ObsData.txt', R_obs,fmt='%4.5f')  
    
    if reduce2d:
        np.savetxt(path  + '/VRTeCoord.txt', 
                    # nodes[:,[0,2]],  
                    np.c_[nodes[imin,0],nodes[imin,2]],  
                    delimiter='\t', fmt='%0.4f')  
    else:
        np.savetxt(path  + '/VRTeCoord.txt', 
                    np.c_[nodes[imin,0],nodes[imin,2],nodes[imin,1]],  
                    delimiter='\t', fmt='%0.4f')  



def create_tmp_MALM_project(k_MALM,i):
    
    if not os.path.exists(k_MALM.dirname + str(i)):
        os.makedirs(k_MALM.dirname + str(i))
        
    k_MALM_tmp = Project(k_MALM.dirname + str(i), typ="R3t")
    k_MALM_tmp = k_MALM
    
    return k_MALM_tmp


def par_fwd_vrte(k_MALM,k_ERT,imin):
    R_icsd = []
    vrtefwd_args = partial(_run_fwd_vrte_par, 
                              k_MALM
                              )            
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results_fwd_vrte = pool.map(vrtefwd_args, 
                                      list(imin)
                                      )
        print(f"x= {imin}, PID = {os.getpid()}")
    print(results_fwd_vrte)
    return results_fwd_vrte
    
    
def _run_fwd_vrte_par(k_MALM,imin):

    # build top part of the matrix A (time consuming part)
    # compute response for each possible source
    # Rsim = np.zeros((k_MALM.sequence.shape[0], len(imin)))
    # k_MALM.mesh.df['res0'] = k_ERT.meshResults[0].df['Resistivity(ohm.m)']
    # k_MALM.mesh.df['res0'] = 100
    # with io.capture_output() as captured:
    # print('\n=====', i, '++++', len(imin))
    # k_MALM.param['node_elec'][1][1] = imin
    k_MALM.forward()
    # Rsim[:,i] = k_MALM.surveys[0].df['resist'].values
    # R_icsd.append(k_MALM.surveys[0].df['resist'].values)
    # R_icsd = k_MALM.surveys[0].df['resist'].values
    R_icsd = imin
    return R_icsd

def fwd_vrte(k_MALM,k_ERT,imin, rho0=None):
    # build top part of the matrix A (time consuming part)
    # compute response for each possible source
    Rsim = np.zeros((k_MALM.sequence.shape[0], len(imin)))
    
    if rho0:
        k_MALM.mesh.df['res0'] = rho0
    else:
        print('use ERT inverted as initial model for Green fcts simulation')
        k_MALM.mesh.df['res0'] = k_ERT.meshResults[0].df['Resistivity(ohm.m)']
        # k_MALM.mesh.df['res0'] = 100
                        
    R_icsd = []
    
    # k_MALM.mesh.node[imin]
    # nodes = k_MALM.mesh.node

    k_ERT.param['node_elec']    
    k_MALM.param['node_elec']    
    # k_MALM.mesh.node[80]
    # k_MALM.param['node_elec'][0] = k_ERT.param['node_elec'][0] 

    # k_MALM.param['node_elec'][1] = k_MALM.param['node_elec'][1]-1

    # d = 80
    for i, d in enumerate(imin):
        # print(d)
        # with io.capture_output() as captured:
        print('\n=====', i, '++++', len(imin))
        # k.param['node_elec'][1][1] = imin[s] # the source to be found (CHANGE HERE)
        k_MALM.surveys[0].df[['a','b','m','n']]
        k_MALM.param['node_elec'][1][71] = d
        # k_MALM.mesh.node[k_MALM.param['node_elec'][1]-1]
        k_MALM.mesh.node[k_MALM.param['node_elec'][1][71]]
        k_MALM.forward()
        # plt.plot(k_MALM.surveys[0].df['resist'].values)
        # plt.show()
        Rsim[:,i] = k_MALM.surveys[0].df['resist'].values
        # plt.plot(Rsim[:,i])
        R_icsd.append(k_MALM.surveys[0].df['resist'].values)
        
        # RobsSi[:,i] = k.surveys[0].df['resist'].values

        
    return R_icsd, Rsim
    
def filter_sequence(k_MALM,**kwargs):
    
    # set default values
    a=72
    b=65
    m=None
    
    if 'a' in kwargs:
        a = kwargs['a']
    if 'b' in kwargs:
        b = kwargs['b']
    if 'm' in kwargs:
        m = kwargs['m']
        
    df_new = k_MALM.surveys[0].df  
    
    if m is not None:
        condition = (df_new['a']!=a) | (df_new['b']!=b) | (df_new['m']!=m)
    else:
        condition = (df_new['a']!=a) | (df_new['b']!=b)

        
    df_new.drop(df_new[condition].index, inplace=True)
    Robs = df_new['resist']
    
    return Robs
    
def filter_sequence_rec(k_MALM):
    
    df_new = k_MALM.surveys[0].df  
    condition = (df_new['a']!=72)
    
    np.sum(condition)
    df_new.drop(df_new[condition].index, inplace=True)
    Robs = df_new['resist']
    
    return Robs

    
def create_MALM_project(path,filename,**kwargs):
    
    # root_path = os.path.split(path)[0] + '/'
    # root_path = str(pathlib.Path(path).parent) +'/'
    # try:
    #     root_path = str(pathlib.Path(path).parent.parent) +'/'
    # except:
    #     pass

    root_path = os.path.join(*path.split('/')[0:2]) +'/'
    
    percent_rec=10
    if 'percent_rec' in kwargs:
        percent_rec = kwargs['percent_rec']
        
    # filter_seq={'bool_filt':False,'a':72,'b':65,'m':71}
    filter_seq=False
    if 'filter_seq' in kwargs:
        filter_seq = kwargs['filter_seq']

    filter_seq_rec=True
    if 'filter_seq_rec' in kwargs:
        filter_seq_rec = kwargs['filter_seq_rec']
        
    if not os.path.exists(path + filename):
        os.makedirs(path + filename)


    k_MALM = Project(path + filename, typ="R3t")
    # k.createSurvey('/home/ben/Documents/Ongoing/Rhizotron_PRD/PRD/rawData/'+ name + '.csv', ftype="Syscal")
    k_MALM.createSurvey(root_path + 'rawData/' + filename + '.csv', ftype="Syscal")
    
    # # bug when importing test sequence replacing 71 instead of 72!!
    # because one missing electrode in the sequence so reduced to 71!
    # if k_MALM.surveys[0].df['a'].all() == '71':
    #     k_MALM.surveys[0].df['a'] = k_MALM.surveys[0].df['a'].replace('71','72')
    k_MALM.surveys[0].df['a'] 
    len(k_MALM.surveys[0].df['a'])
    
    k_MALM.importElec(root_path + "mesh/elecsXYZ.csv")
    
    # k_MALM.sequence
    # k_MALM.elec
    # print(len(k_MALM.surveys[0].df['a']))
    
    # k_MALM.importMesh(root_path + "mesh/mesh_rhizo_resipy.msh")
    k_MALM.importMesh(root_path + MESH,
                        # elec=k_MALM.elec[['x','y','z']].to_numpy(),
                        # node_pos=np.arange(72)
                      )
    
    # k_MALM.importElec(root_path + "mesh/elecsXYZ.csv")

    nodes = k_MALM.mesh.node
    
    k_MALM.elec
    # dist = cdist(grid, nodes[:])
    # imin = np.argmin(dist, axis=1)

    
    # newnode_elec = []  
    # for i in k_MALM.param['node_elec'][0]:
    #     newnode_elec.append(str(int(i)-1))
    
    
    # k_MALM.param['node_elec'][0] =  np.array(newnode_elec)
    # k_MALM.param['node_elec'][1] =  k_MALM.param['node_elec'][1] - 1


    # elecsMALM_nodes = k_MALM.param['node_elec'][1]
    # nodes[71]
    
    # len(nodes[elecsMALM_nodes])
    # k_MALM.surveys[0].df[['a','b','m','n']]
    # test = k_MALM.surveys[0].df[['resist']]


    try:
        k_MALM.filterRecip(percent=percent_rec)
        # k_MALM.filterRecip(percent=10)
    except:
        print('can''t find reciprocals')
        # pass
    k_MALM.surveys[0].df[['a','b','m','n']]= k_MALM.surveys[0].df[['a','b','m','n']].astype('int64')

       
    # test = k_MALM.surveys[0].df[['resist']]

    if (filter_seq & filter_seq_rec):
        raise ValueError('impossible combinaison of filters')
        
    if filter_seq:
        Robs = filter_sequence(k_MALM,**kwargs)
    elif filter_seq_rec:
        Robs = filter_sequence_rec(k_MALM)
    else:
        Robs =  k_MALM.surveys[0].df['resist']
        # Robs =  k_MALM.surveys[0].df['rho']

        
    # plt.scatter(k_MALM.surveys[0].df['resist'].index, k_MALM.surveys[0].df['resist'])
    # plt.scatter(k_MALM.surveys[0].df['Rho'].index, k_MALM.surveys[0].df['Rho'])


    k_MALM.surveys[0].df[['a','b','m','n']]= k_MALM.surveys[0].df[['a','b','m','n']].astype('str')


    # k_MALM.surveys[0].df['a'] 
    # k_MALM.surveys[0].df['b'] 
    len(k_MALM.surveys[0].df['a'])
    # k.showErrorDist()
    # k.filterManual(index=0, attr="resist" # interactive plot)
    
    if 'ax' in kwargs:
        Robs.plot(ylabel=r'R ($\Omega$)', xlabel='n #', 
                  legend='A=(stem of soil), B=71, M=64',
                  style='o-', label=filename, ax=kwargs['ax'])
    Robs.rename(filename, inplace=True)
    len(Robs)
    return k_MALM, Robs


def create_source_grid(mesh, k, nVRTe_rows=5, nVRTe_cols=4, show=False, **kwargs):
    
    
    # define rhizo size
    
    # width = (490 - 10 -10)*1e-3
    # height = (500)*1e-3

    # # dx = (1/nVRTe_rows)/2.2
    # offset = 0.0075
    # dx = (0.47-offset*2)/(nVRTe_rows-1)
    # dz = dx + 0.002
    
    #create a regular grid of possible sources
    # gridx, gridz = np.meshgrid(min(mesh['X'])+offset + np.arange(nVRTe_rows)*dx, 
    #                             min(mesh['Z'])+offset + np.arange(nVRTe_cols)*dz
    #                             )
    
    
    # minx_elec = 0.06
    # elec_spac = 0.05
    # maxx_elec = 0.41
    # miny_elec = 0.085
    # maxy_elec = 0.487
    # offset = elec_spac/2
    # nVRTe_rows = 7
    
    # gridx, gridz = np.meshgrid(minx_elec+offset + np.arange(nVRTe_rows)*elec_spac, 
    #                            miny_elec+offset + np.arange(nVRTe_rows+1)*elec_spac
    #                             )
    
    
    # minx_elec = 0.06
    elec_spac = 0.05
    # maxx_elec = 0.41
    miny_elec = 0.085
    maxy_elec = 0.487
    offset = elec_spac/2
    # nVRTe_rows = 7
    eps=1e-4
    
    gridx, gridz = np.meshgrid(
                                np.linspace(0+eps,47*1e-2-eps,nVRTe_cols), 
                                np.linspace(0+eps,50*1e-2-eps,nVRTe_rows), 
                               # np.arange(nVRTe_cols)*elec_spac*2 + offset
                                )
    
    # (0.45-offset*2)/10
    # gridz = gridz*-1
    gridx = gridx.flatten()
    
    if gridx.min()<min(mesh['X']):
        print('error min X')
    if gridx.max()>max(mesh['X']):
        print('error max X')
    if gridz.min()<min(mesh['Z']):
        print('error min Z')
    if gridz.max()>max(mesh['Z']):
        print('error max Z')
        
    # mesh['Y'].max()
    gridz = gridz.flatten()
    grid = np.c_[gridx, 
                 np.zeros(len(gridx))*mesh['Y'].max(),
                 gridz]
    nodes = k.mesh.node
    
    # snap grid points to mesh nodes
    dist = cdist(grid, nodes[:])
    imin = np.argmin(dist, axis=1)


    elecsMALM_nodes = k.param['node_elec'][1]
    k.elec 
        
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.scatter(nodes[imin,0],nodes[imin,2])
    # ax.scatter(gridx,gridz)
    ax.set_xlim([0,0.47])
    ax.set_ylim([0,0.5])
    
    # max(gridx)
    # min(gridx)
    # len(imin)
    # check that nodes does not belong to electrodes nodes
    elecsMALM = k.surveys[0].elec
    # plt.scatter(nodes[elecsMALM_nodes,0],
    #             nodes[elecsMALM_nodes,2])
    # plt.scatter(elecsMALM['x'],elecsMALM['z'])    
    
    
    interection_danger = np.intersect1d(elecsMALM_nodes,imin)
    
    if len(interection_danger)>1:
        # raise ValueError('Intersection not null')
        print('ValueError Intersection not null')
        
    # def check_VRTenodesVS_elecsnode(nodes,imin,elecsMALM):
        
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
        # print(i)
        if i==(90+72):
            stri += str(i+1)
        else:
            stri += str(i+1) + ','
    n = text_file2.write(stri)
    text_file2.close()

        
        
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.scatter(nodes[imin,0],nodes[imin,2],nodes[imin,1])
    
    # ax.scatter(nodes[imin,0],nodes[imin,2],nodes[imin,1])

    if show==True:
        ax = pv.Plotter()
        # k.showMesh(ax=ax,pvshow=False,
        #                 xlim=[mesh['X'].min(),mesh['X'].max()],
        #                 ylim=[mesh['Y'].min(),mesh['Y'].max()], 
        #                 )
        ax.add_points(nodes[imin,:])
        ax.show(auto_close=True)
    
    return imin, nodes, grid

def set_MALM_root_seq(k):
    k.sequence = k.surveys[0].df[['a','b','m','n']]
    return k

def _run_invert_ERT_TL(path, filepaths, folderName, regType=1,
                  recip=10,
                  **kwargs):
    
    k = Project(path + folderName, typ="R3t")
    # k.createSurvey('/home/ben/Documents/Ongoing/Rhizotron_PRD/PRD/rawData/'+ name + '.csv', ftype="Syscal")
    k.createTimeLapseSurvey(filepaths, ftype="Syscal")
    # k.showError(index=0)
    k.filterRecip(percent=recip)
    k.filterRecip(index=-1, percent=recip)

    # k.showErrorDist()
    # k.filterManual(index=0, attr="resist" # interactive plot)
    #k.showPseudo()
    k.importElec(path + "mesh/elecsXYZ.csv")
    #k.setElec(elec)
    # k.filterElec([28])
    # k.showPseudo()
    # k.filterManual(index=0, attr="resist" # interactive plot)
    k.importMesh(path + MESH)
    # k.showMesh()
    # k.setStartingRes(regionValues={1: 100.0}, 
    #                  zoneValues={1: 1}, 
    #                  fixedValues={1: False}, 
    #                  ipValues={1: 0.0}) # define region manually using k.addRegion()
    
    k.param['reg_mode'] = regType
    
    if regType==1:
        k.param['num_xy_poly'] = 0 # tells R3t to no crop the mesh
        k.param['zmin'] = -np.inf
        k.param['zmax'] = np.inf
        # k.param['ymin'] = -np.inf
        # k.param['ymax'] = np.inf
        
    # k.err = True # if we want to use the error from the error models fitted before
    k.param['a_wgt'] = 0.01
    k.param['b_wgt'] = recip/100

    k.invert(modErr=False, parallel=False, modelDOI=False)
    k.saveProject(path + folderName + 'backup')   
    return k

    
def invert_ERT_TL(path='.',
                  files='',
                  regType=1,
                  recip=10,
                  reprocessed=False,
                  **kwargs):
        
    folder_ERT_TL = 'inversion_ERT_TL/'
    if 'idfileNames' in kwargs:
        idTL = kwargs['idfileNames']
    fileNames = str([*idTL]).replace(", ", "_")
    fileNames = fileNames.replace("[", "")
    fileNames = fileNames.replace("]", "")

    fileNames += 'Reg' + str(regType) 
    
    if not os.path.exists(path + folder_ERT_TL + fileNames):
        os.makedirs(path + folder_ERT_TL + fileNames)
    
    filepaths = []
    for f in files:
        filepaths.append(path + 'rawData/' + f + '.csv')
       
    if os.path.exists(os.path.join(path + folder_ERT_TL + fileNames + "backup.resipy")):
        if reprocessed:
            k = _run_invert_ERT_TL(path, filepaths, folder_ERT_TL + fileNames, regType=regType,
                              recip=recip,
                              **kwargs)
        else:
            print('reload file from *.resipy backup')
            k = Project(dirname=path+folder_ERT_TL + fileNames, typ='R3t')
            k.loadProject(os.path.join(path,folder_ERT_TL + fileNames  + "backup.resipy"))
    else:
        k = _run_invert_ERT_TL(path, filepaths, folder_ERT_TL + fileNames, regType=regType,
                          recip=recip,
                          **kwargs)
            
    return [k]


def plot_scatter_MALM(imaging_path,df_MALM,k,
                        ERT_log, f_MALM,
                        **kwargs):
        ''' Loop over MALM files and plot '''
        # df = pd.concat(R_obs_stck, axis=1, keys=[s.name for s in R_obs_stck])
        
        elecsMALM = k[0].surveys[0].elec
        elecsMALM = elecsMALM.drop(index=[64,70,71], axis=0)
        
        xelecs = elecsMALM['x'] 
        yelecs = elecsMALM['z']
                
        for i, f in enumerate(df_MALM.columns[:-1]):
           
            fig, ax = plt.subplots(
                                   constrained_layout=False)
    
            # lgd += 'Stem-Stem (Ohm)'
            if 'log_scale' in kwargs:
                cmap = ax.scatter(x=xelecs,y=yelecs,c=df_MALM[f[0]].to_numpy(), s=200,
                                  norm=mpl.colors.LogNorm(),
                                  **kwargs)
                cbar = plt.colorbar(cmap, ax=ax)
                cbar.set_label(r'R ($\Omega$)')#, rotation=270)
            else:
                cmap = ax.scatter(x=xelecs,y=yelecs,c=df_MALM[f[0]].to_numpy(), s=200,
                                  **kwargs)
                cbar = plt.colorbar(cmap, ax=ax)
                cbar.set_label(r'R ($\Omega$)') #, rotation=270)
                
            # Loop for annotation of all points
            for i_label,j_label in enumerate(elecsMALM['label']):
                plt.annotate(j_label,(xelecs.iloc[i_label],yelecs.iloc[i_label]))
    
            ax.set_xlabel('x (m)')
            ax.set_ylabel('z (m)')
    
            ax.set_title(str(f[0]))
            
            # with open(imaging_path + str(f[0]) + '.obj', 'wb') as file:
            #     pickle.dump((fig, ax), file)
            plt.savefig(imaging_path + f_MALM[i] + '/' + str(f[0]) + '.png',
                        dpi=450)
            
        return df_MALM
    
def plot_scatter_MALM_diff(imaging_path,df_MALM,k_indiv_merged,
                           selected_files_MALM,
                           ERT_log, background_diff=True,
                           **kwargs):
        ''' Loop over MALM files and compute pair differences '''
        # df = pd.concat(R_obs_stck, axis=1, keys=[s.name for s in R_obs_stck])
        
        # selected_files_MALM = df_MALM.columns
        
        elecsMALM = k_indiv_merged[0].surveys[0].elec
        elecsMALM = elecsMALM.drop(index=[64,70,71], axis=0)
        
        xelecs = elecsMALM['x'] 
        yelecs = elecsMALM['z']
        
        df_MALM['diff_level'] = False
        df_MALM = df_MALM.set_index('diff_level', append=True).unstack('diff_level')

        for injType in ['Stem', 'Soil']:
            cond_same_nb_returns= ~ERT_log['Name'].str.contains('8returns')
            if injType == 'Stem':
                subkey_name = 'diff_stem'
                selected_files_MALM_injType = ERT_log['Name'][(ERT_log['injection']=='Stem')&
                                                              (cond_same_nb_returns)
                                                              ].to_list()
            else:
                subkey_name = 'diff_soil'
                # cond_same_nb_returns= ~ERT_log['Name'].str.contains('8returns')
                selected_files_MALM_injType = ERT_log['Name'][(ERT_log['injection']=='Soil') & 
                                                              (ERT_log['method']=='MALM')  &
                                                              (cond_same_nb_returns)
                                                              ].to_list()

            
            # df_diff = pd.DataFrame()
            for i, f in enumerate(selected_files_MALM_injType):
                if (i < len(selected_files_MALM_injType)-1):
                    fig, ax = plt.subplots(
                                           constrained_layout=False)
                    
                    d1=ERT_log[ERT_log['Name'] == selected_files_MALM_injType[i+1]]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
        
                    if background_diff:
                        id_file = 0
                        d0=ERT_log[ERT_log['Name'] == selected_files_MALM_injType[0]]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
                    else:
                        id_file = i
                        d0=ERT_log[ERT_log['Name'] == selected_files_MALM_injType[0]]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
        
        
                    diff_keyname = (d0 + ' - ' + d1)
                    
                    # df_MALM[diff_keyname] = np.NaN
                    
                    d_back=ERT_log[ERT_log['Name'] == selected_files_MALM_injType[id_file]]['datetime'].to_list()[0]
                    d_forth=ERT_log[ERT_log['Name'] == selected_files_MALM_injType[i+1]]['datetime'].to_list()[0]
                    key_datetime =  d_forth
                    
                    df_MALM[(key_datetime,diff_keyname,subkey_name)] = np.NaN
                    df_MALM[(key_datetime,diff_keyname,subkey_name)] = (df_MALM[(d_back,injType,False)] - 
                                                                        df_MALM[(d_forth,injType,False)]
                                                                        ).to_numpy()
        
                        
                    # lgd += 'Stem-Stem (Ohm)'
                    if 'log_scale' in kwargs:
                        cmap = ax.scatter(x=xelecs,y=yelecs,c=(df_MALM[(key_datetime, 
                                                                           diff_keyname,
                                                                           subkey_name)
                                                                          ]), s=200,
                                          norm=mpl.colors.LogNorm(),
                                          **kwargs)
                        cbar = plt.colorbar(cmap, ax=ax)
                        cbar.set_label(subkey_name + '(Ohm)')#, rotation=270)
                    else:
                        cmap = ax.scatter(x=xelecs,y=yelecs,c=(df_MALM[(key_datetime, 
                                                                           diff_keyname,
                                                                           subkey_name)
                                                                          ]), s=200,
                                          **kwargs)
                        cbar = plt.colorbar(cmap, ax=ax)
                        cbar.set_label(subkey_name+ '(Ohm)') #, rotation=270)
                    ax.set_xlabel('x (m)')
                    ax.set_ylabel('z (m)')
        
                    ax.set_title((diff_keyname,subkey_name))
                    
                    plt.savefig(imaging_path + f + '/' + f + '_diff.png',
                                dpi=450)
                    # with open(imaging_path + f + '_diff' + '.obj', 'wb') as file:
                    #     pickle.dump((fig, ax), file)
    
        return df_MALM
                

def plot_scatter_MALM_diff_stem_soil(imaging_path,df_MALM,k_indiv_merged,
                                     selected_files_MALM,
                                     ERT_log, **kwargs):

    # df = pd.concat(R_obs_stck, axis=1, keys=[s.name for s in R_obs_stck])
    
    elecsMALM = k_indiv_merged[0].surveys[0].elec
    elecsMALM = elecsMALM.drop(index=[64,70,71], axis=0)
    
    xelecs = elecsMALM['x'] 
    yelecs = elecsMALM['z']
            
    # selected_files_MALM = df_MALM.columns[0]
    
    for i, f in enumerate(selected_files_MALM):
        j = int(math.floor(i/2))    
        if (i % 2) == 0:
            fig, ax = plt.subplots(
                                   constrained_layout=False)
            
            df_MALM['diff_stem_soil' + str(j)] = np.NaN
            
            d0=ERT_log[ERT_log['Name'] == f]['datetime'].to_list()[0]
            d1=ERT_log[ERT_log['Name'] == selected_files_MALM[i+1]]['datetime'].to_list()[0]

            df_MALM['diff' + str(j)] = d0 - d1
            df_MALM['diff' + str(j)].isnull().values.any()
            # df_MALM[f].isnull().values.any()
            
            d0=ERT_log[ERT_log['Name'] == f]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
            d1=ERT_log[ERT_log['Name'] == selected_files_MALM[i+1]]['datetime'].to_list()[0].strftime('%b %d,%H:%M')
                        
            lgd = (d0 + ' - ' + d1)

            if 'log_scale' in kwargs:
                cmap = ax.scatter(x=xelecs,y=yelecs,c=(df_MALM['diff' + str(j)]), s=200,
                                  norm=mpl.colors.LogNorm()
                                  **kwargs)
                cbar = plt.colorbar(cmap, ax=ax)
                cbar.set_label('|Stem-Soil| (Ohm)')#, rotation=270)
            else:
                cmap = ax.scatter(x=xelecs,y=yelecs,c=(df_MALM['diff' + str(j)]), s=200,
                                  **kwargs)
                cbar = plt.colorbar(cmap, ax=ax,
                                    )
                cbar.set_label('Stem-Soil (Ohm)')#, rotation=270)
                
            ax.set_xlabel('x (m)')
            ax.set_ylabel('z (m)')

            ax.set_title(lgd)
            lgd += 'Stem-Soil (Ohm)'
            
            plt.savefig(imaging_path + f +'/' + f + '_SS_diff.png',
                        dpi=450)
            # with open(imaging_path + f + '_SS_diff' + '.obj', 'wb') as file:
            #     pickle.dump((fig, ax), file)
                
            # plt.savefig(imaging_path + str(f[0]) + '_SS_diff.png',
            #             dpi=450)
            # with open(imaging_path + str(f[0]) + '_SS_diff' + '.obj', 'wb') as file:
            #     pickle.dump((fig, ax), file)
                
    return df_MALM

            
            
def plot_MALM_deprecated(k_MALM,nodes,imin,sol,show=True):
    ax = pv.Plotter(notebook=False,window_size=([2048, 1536]))
    k_MALM.showMesh(ax=ax,pvshow=False,
                    xlim=[k_MALM.mesh.df['X'].min(),k_MALM.mesh.df['X'].max()],
                    ylim=[k_MALM.mesh.df['Y'].min(),k_MALM.mesh.df['Y'].max()], 
                    )
    
    points =  np.c_[nodes[imin,0],-np.ones(len(nodes[imin,2]))*0.005,nodes[imin,2]]
    cloud = pv.wrap(points)
    ax.add_points(
                    cloud,
                    scalars=sol.x,
                    render_points_as_spheres=True,
                    point_size=25,
                    show_scalar_bar=False
                    # opacity=sol.x*10,
                )
    _ = ax.add_scalar_bar('csd', interactive=False, vertical=True,
                            title_font_size=25,
                            label_font_size=20,
                            outline=True, fmt='%10.5f')
    ax.save_graphic(k_MALM.dirname + "MALM.svg")
    ax.save_graphic(k_MALM.dirname + "MALM.eps")
    
    if show:
        ax.show()
    else:
        return ax
    
def to_mpl(mesh):
    
    z = 7.500e-02
    x = None
    y = None
    slices = mesh.ctp().slice_orthogonal(x=x, y=y, z=z, 
                               generate_triangles=True, 
                               contour=True)
    
    
    cmap = plt.cm.get_cmap("viridis", 4)
    slices = mesh.slice_orthogonal()
    slices.plot(cmap=cmap)

    for slc in slices:
        pts = slc.points
        tri = slc.faces.reshape((-1,4))[:, 1:]
        val = slc.active_scalar
        
    data = mesh.ctp().slice("x", generate_triangles=True)

    # Grab data for Matplotlib
    x = data.points
    tri = data.faces.reshape((-1,4))[:, 1:]
    u = data.active_scalar
    
    # Plot it yo!
    plt.tricontourf(x[:,0], x[:,1], tri, u)
    plt.gca().set_aspect('equal')


def plot_ERT(k, vmin=0, vmax=300, attr="Resistivity(log10)", index=0, path='',
             ext='.png', ax=None,
             show=True, showElec=True,
             **kwargs):
    '''
    '''
    #%%
    import vtk
    color_map = 'jet'
    titlecbar = r'ER' + '\n' + r'($\Omega$.m)'
    if "Resistivity(log10)" in attr:
        titlecbar = r'log10(ER)' + '\n' + r'($\Omega$.m)'

    if "difference(percent)" in attr:
        color_map = 'bwr'
        titlecbar = r'$\Delta ER$' + '\n' + r'$\%$'
    
    fi = ''
    if 'fi' in kwargs:
        fi = kwargs['fi']
    mesh = k.mesh.df

    # pl.colorbar_orientation = 'vertical'  
    # pl.font.size = 20  
    # pl.font.label_size = 10  
    
    if ax is None:
        ax = pv.Plotter(notebook=True,window_size=([2048, 1536])
                        )
        # ax = pv.Plotter(notebook=True,window_size=([1024, int(1536/2)])
        #                 )
        
   
    sargs = dict(height=0.35, 
                 width=0.05,
                 vertical=True, 
                 position_x=0.035, 
                 position_y=0.6,
                 title_font_size=30,
                 label_font_size=30,
                 title=titlecbar,
                 )

    k.showResults(
                  index=index,
                  attr=attr,
                  ax=ax, 
                  vmin=vmin, vmax=vmax,
                  color_map=color_map,
                  background_color='white',
                  pvgrid = False,
                  use_pyvista=True,
                  pvshow=False,
                  xlim=[mesh['X'].min(),mesh['X'].max()],
                  ylim=[mesh['Y'].min(),mesh['Y'].max()], 
                  zlim=[mesh['Z'].min(),mesh['Z'].max()],  
                   scalar_bar_args=sargs,
                  pvdelaunay3d=True,
                  # color_bar=True,
                  # hor_cbar=False,
                  # )
                   **kwargs)
    
    
    # ax.camera_position = 'xz'
    # ax.camera.azimuth = 15
    # ax.camera.elevation = 15
    # # ax.camera.roll += 10



    # pl.show_bounds(all_edges=True)
    if showElec:
        ax.view_xz(negative=True)
    else:
        ax.view_xz()
    
    ax.set_scale(1.1, 1.1, 1.1)
    ax.show_bounds(grid=False, location='origin', 
                   all_edges=False, 
                   font_size=30,
                   xlabel='X (m)',
                   zlabel='Z (m)',
                   show_yaxis=False,
                   # ticks='outside',
                   use_2d=False)


    # pl.scalar_bars
    # pl.view_xz()
    #actor = pl.show_grid()
    # actor = pl.show_bounds(grid='back', location='front',
    #                         all_edges=True, color='black',use_2d=True)
    # pl.remove_scalar_bar()
    # pl.add_scalar_bar()
    # pl.scalar_bars.vertical = True
    
    # pv.global_theme.colorbar_orientation = 'vertical'  
    # pv.global_theme.font.size = 20  
    # pv.global_theme.font.label_size = 10  
    
    # ax.colorbar_orientation = 'vertical'  
    # pl.SetUse2DMode(True)
    # pl.set_scale(1.1, 1.1, 1.1)

    # pl.font.size = 20  
    # pl.font.label_size = 10  
    
    
    
    if 'png' in ext:
        
        hsize = 1150
        # pl.ren_win.OffScreenRenderingOn()
        # pl.enable_anti_aliasing()
        # pl.screenshot('../src/image/paper3d/mesh-types.jpg', transparent_background=True,
        #               window_size=[hsize, int(hsize*pl.window_size[1]/pl.window_size[0])])
        # pl.ren_win.SetSize([1000, 800])
        # pl.ren_win.OffScreenRenderingOff()
        # pl.ren_win.Render()

        ax.screenshot(k.dirname + '/' + attr + '_time' + str(index) + fi + ".png",
                      window_size=[hsize, hsize])
    if 'svg' in ext:
        ax.save_graphic(k.dirname + '/' + attr + '_time' + str(index) + fi + ".svg")  
        # pl.save_graphic(k.dirname + '/' + attr + '_t' + str(index) + ".eps")  

    #%%
    # pl.window_size = window_size
    
    # k.meshResults[0].
    
    
    # to_mpl(pl.mesh)
    

    # if show:
        # ax = pv.Plotter(notebook=True,window_size=([2048, 1536])
        #                 )
        
    # k.showResults(
    #               index=index,
    #               attr=attr,
    #               ax=ax, 
    #               vmin=vmin, vmax=vmax,
    #               color_map=color_map,
    #               background_color='white',
    #               pvgrid = True,
    #               use_pyvista=True,
    #               pvshow=False,
    #               xlim=[mesh['X'].min(),mesh['X'].max()],
    #               ylim=[mesh['Y'].min(),mesh['Y'].max()], 
    #               zlim=[mesh['Z'].min(),mesh['Z'].max()], 
    #               **kwargs)
    
    
    # # ax.camera_position = 'yz'
    # ax.camera.azimuth = -125
    # ax.camera.elevation = -10
    # # ax.camera.roll += 10


    # pl.show_bounds(all_edges=True)
    # pl.view_xz(negative=True)
    # pl.scalar_bars
    #ax.view_xz()
    #actor = pl.show_grid()
    # actor = pl.show_bounds(grid='back', location='front',
    #                         all_edges=True, color='black',use_2d=True)
    # pl.remove_scalar_bar()
    # pl.add_scalar_bar()
    # pl.scalar_bars.vertical = True
    
    # pv.global_theme.colorbar_orientation = 'vertical'  
    # pv.global_theme.font.size = 20  
    # pv.global_theme.font.label_size = 10  
    
    # ax.colorbar_orientation = 'vertical'  
    # pl.SetUse2DMode(True)
    # pl.set_scale(1.1, 1.1, 1.1)

    # pl.font.size = 20  
    # pl.font.label_size = 10  
    
    
    
    # if 'png' in ext:
        
    #     # hsize = 2000
    #     # pl.ren_win.OffScreenRenderingOn()
    #     # pl.enable_anti_aliasing()
    #     # pl.screenshot('../src/image/paper3d/mesh-types.jpg', transparent_background=True,
    #     #               window_size=[hsize, int(hsize*pl.window_size[1]/pl.window_size[0])])
    #     # pl.ren_win.SetSize([1000, 800])
    #     # pl.ren_win.OffScreenRenderingOff()
    #     # pl.ren_win.Render()

    #     ax.screenshot(k.dirname + '/' + attr + '_t' + str(index) + ".png")
    # if 'svg' in ext:
    #     ax.save_graphic(k.dirname + '/' + attr + '_t' + str(index) + ".svg")  
    #     # pl.save_graphic(k.dirname + '/' + attr + '_t' + str(index) + ".eps")  

    # pl.window_size = window_size
    
    # k.meshResults[0].
    
    
    # to_mpl(pl.mesh)
    

    # if show:
    #     # ax.show()
    # else:
    return ax

def invert_ERT(path='.',
               filename='',
               recip=10,
               **kwargs):

    
    if not os.path.exists(path + 'inversionERT/' + filename):
        os.makedirs(path + 'inversionERT/' + filename)


    k = Project(path + 'inversionERT/' + filename, typ="R3t")
    # k.createSurvey('/home/ben/Documents/Ongoing/Rhizotron_PRD/PRD/rawData/'+ name + '.csv', ftype="Syscal")
    k.createSurvey(path + 'rawData/' + filename + '.csv', ftype="Syscal")
    # k.showError(index=0)
    k.filterRecip(percent=recip)
    # k.showErrorDist()
    # k.filterManual(index=0, attr="resist" # interactive plot)
    # k.showPseudo()
    k.importElec(path + "mesh/elecsXYZ.csv")
    # k.filterElec([28])
    # k.filterManual(index=0, attr="resist" # interactive plot)
    # k.importMesh(path + "mesh/mesh_rhizo_resipy.msh")
    k.importMesh(path + MESH)
    # k.createMesh()
    # k.showMesh()
    k.setStartingRes(regionValues={1: 100.0}, 
                     zoneValues={1: 1}, 
                     fixedValues={1: False}, 
                     ipValues={1: 0.0}) # define region manually using k.addRegion()
    k.param['num_xy_poly'] = 0 # tells R3t to no crop the mesh
    k.param['zmin'] = -np.inf
    k.param['zmax'] = np.inf

    # k.err = True # if we want to use the error from the error models fitted before
    k.param['a_wgt'] = 0.001
    k.param['b_wgt'] = recip/100
    # k.saveSequence('ERTseq.csv')
    k.surveys[0].df[['a','b','m','n']].to_csv('seqERT.csv', index=False)

    k.invert(modErr=False, parallel=False, modelDOI=False)
    k.saveProject(path + 'inversionERT/' + filename + 'bkp')

    return k


def getR2out(imaging_path, files, typ='R2'):
    getR3out(imaging_path, files, typ='R2')
    
def getR3out(imaging_path, files, typ='R3t'):
    """Reat the .out file and parse its content.
    
    Returns
    -------
    Dataframe with the dataset name, and the RMS decrease for each iteration.
    """
    
    DF_RMS = []
    for i, file_survey in enumerate(files):
        fname = os.path.join(imaging_path, 'inversionERT') + '/' + file_survey + '/invdir/' + typ + '.out'
        
        with open(fname, 'r') as f:
            lines = f.readlines()
        name = ''
        idataset = 0
        iiter = 0
        resRMS = np.nan
        phaseRMS = np.nan
        read = np.nan
        rejected = np.nan
        irow = 0
        df = pd.DataFrame(columns=['name', 'dataset', 'iteration', 'resRMS',
                                   'phaseRMS', 'read', 'rejected', 'success'])
        for x in lines:
            success = 'N/A'
            line = x.split()
            if len(line) > 1:
                if line[0] == 'Iteration':
                    iiter += 1
                elif ('measurements' in x) & ('read:' in x):
                    print(line)
                    read = int(line[5])
                    # rejected = int(line[5])
                elif line[0] == 'Final':
                    resRMS = float(line[3])
                    df.loc[irow, :] = [name, idataset, iiter, resRMS, phaseRMS,
                                       read, rejected, success]
                    irow += 1
                elif line[0] == 'FATAL:':
                    resRMS = np.nan
                elif line[0] == 'Processing':
                    iiter = 0
                    idataset += 1
                    # if idataset <= len(surveys[i]):
                    #     name = i
                    # else:
                    # name = 'dataset{:03.0f}'.format(idataset)
                    name = file_survey
        df = df.apply(pd.to_numeric, errors='ignore').reset_index(drop=True)
        DF_RMS.append(df)

    DF_RMS = pd.concat(DF_RMS)
        
    return DF_RMS
    
    