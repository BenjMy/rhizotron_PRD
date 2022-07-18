#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:19:06 2022

@author: ben
"""

#!/usr/bin/env python
# coding: utf-8

# In[1]:
import matplotlib.pyplot as plt
from pyPRD import processing as proc
import numpy as np
import pyvista as pv
from resipy import Project
import os

imaging_path = '../../imaging_ERT_MALM/'

# filenames_ERT = ['PRD_ERT_1_130522', 'PRD_ERT_1_170522']
# filenames_MALM = ['PRD_MALM_1_130522', 'PRD_MALM_1_170522']

filenames_ERT = ['PRD_ERT_test0']
filenames_MALM = ['PRD_MALM_test0']
# filenames_MALM = ['PRD_MALM_1_130522_only72_short_TEST']

# In[2]:

k_indiv = []

k = Project(dirname=imaging_path+filenames_ERT[0], typ='R3t')
k.loadProject(os.path.join(imaging_path,
                            'inversionERT',
                            filenames_ERT[0] + "bkp.resipy"
                            )
              )
k_indiv.append(k)


# k_indiv = []
# for i, f in enumerate(filenames_ERT):
#     k_indiv.append(proc.invert_ERT(imaging_path,
#                                    filename=f,
#                                    recip=10)
#                    )
#     proc.plot_ERT(k_indiv[i], vmin=0, vmax=50,
#                   attr="Resistivity(ohm.m)",
#                   index=0)
# proc.plot_ERT(k_indiv[0], vmin=0, vmax=3,
#               attr="Resistivity(log10)",
#               index=0)


# In[ ]:

k_MALM = []
imin_S = [1]  # 0.25,-5e-3,0.25 66
# imin_S = [4,10,30]  # 0.25,-5e-3,0.25 66
# imin_S = [6]  # 0.25,-5e-3,0.25 66

wS = list(np.ones(len(imin_S))/len(imin_S))

for i, f in enumerate(filenames_MALM):
    outMALM = proc.prepare_MALM_synth(imaging_path,
                                      f,
                                      k_indiv[i],
                                      idS=imin_S,
                                      wS=wS,
                                      reduce2d=True,
                                      # filter_seq_rec=True,
                                      filter_seqc=True,
                                      nVRTe_rows=9, nVRTe_cols=9,
                                      )
    k_MALM.append(outMALM[0])
nodes = outMALM[3]
imin = outMALM[1]
R_obs = outMALM[2]
R_icsd = outMALM[4]
R_sim = outMALM[5]

plt.plot(R_obs)
# plt.plot(R_icsd[imin_S[0]])
plt.plot(R_sim[imin_S[0]])

# R_sim= R_sim.T
# np.shape(R_sim[imin_S[0]])
# plt.plot(R_icsd[10])
# points = np.c_[nodes[imin_S, 0], -
#                 np.ones(len(nodes[imin_S, 2]))*0.005, nodes[imin_S, 2]]
# print(points)
# print(nodes[imin[18]])
print(nodes[imin[imin_S]])

# from scipy.spatial.distance import cdist
# tosearch = [[0.3, 0., 0.3]]
# dist = cdist(nodes[imin], tosearch)
# imin_sol = np.argmin(dist, axis=0)
# nodes[imin[27]]
# nodes[imin[66]]

#%%

# from icsd3d.icsd3d_class import iCSD3d as i3d 
# from icsd3d.plotters import mpl_plot
# import os
# os.getcwd()
# path2files= imaging_path + 'inversionMALM/' + f
# icsd=i3d(dirName=path2files +'/')   
# icsd.regMesh='strc'
# icsd.type='2d'
# icsd.obs_err='sqrt' # sqrt choose between constant weight and w = 1/sqrt(abs(obs))
# icsd.wr=1 #weight regularization
# icsd.alphaSxy=False
# icsd.x0_prior=False
# icsd.x0_ini_guess=False # initial guess
# icsd.method_m0='F1'
# icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')
# m0 = icsd.estimateM0(method_m0='F1',show=True)
# # m0 = icsd.estimateM0(method_m0='Pearson',show=True)

# # icsd.clim=[0,0.1]
# icsd.alphax0=1
# sol= icsd.invert(wr=1e-2, pareto=False)

# fig, ax = plt.subplots()
# icsd.showResults(ax=ax)
# ax.scatter(nodes[imin[imin_S],0],nodes[imin[imin_S],2],
#             color='r',marker='v')
# plt.show()


# %%

sol = []
m0 = []
for i, f in enumerate(filenames_MALM):
    sol.append(
                proc.invert_MALM(imaging_path + 'inversionMALM/' + f,
                                 wr=10,
                                 obs_err='sqrt',
                                 x0_ini_guess=False,
                                 alphaSxy=False,
                                 x0_prior=False,
                                 method_m0='F1',
                                 typ='2d',
                                 fname_sim='VRTeSim_icsd.txt',
                                 pareto=False
                                 
                                )
               )

    m0.append(proc.m0_MALM(imaging_path + 'inversionMALM/' + f,
                           method_m0='F1',
                           fname_sim='VRTeSim_icsd.txt',
                           )
              )
    print(nodes[imin[imin_S]])
    proc.plot_m0_MALM(k_MALM[i],m0[i])
    proc.plot_MALM(k_MALM[i], sol[i],show=False,ext=['png'])


# %%

sol = []
m0 = []
for i, f in enumerate(filenames_MALM):
    sol.append(
                proc.invert_MALM(imaging_path + 'inversionMALM/' + f,
                                 wr=1,
                                 obs_err='sqrt',
                                 x0_ini_guess=False,
                                 alphaSxy=False,
                                 x0_prior=False,
                                 method_m0='F1',
                                 typ='2d',
                                 fname_sim='VRTeSim.txt',
                                )
               )

    m0.append(proc.m0_MALM(imaging_path + 'inversionMALM/' + f,
                           method_m0='F1',
                           fname_sim='VRTeSim.txt',
                           )
              )
    print(nodes[imin[imin_S]])
    proc.plot_m0_MALM(k_MALM[i],m0[i])
    proc.plot_MALM(k_MALM[i], sol[i],show=False,ext=['png'])


# %%


# for i, f in enumerate(filenames_MALM):
#     sol = []
#     m0 = []
#     for j, wri in enumerate(np.linspace(1e-2,100,10)):
#         sol.append(proc.invert_MALM(imaging_path + 'inversionMALM/' + f,
#                                     wr=wri)
#                     )
#         m0.append(proc.m0_MALM(imaging_path + 'inversionMALM/' + f,
#                                 method_m0='Pearson'))
#         print(nodes[imin_S])
#         proc.plot_MALM(k_MALM[i],nodes,imin,sol[j],show=True)
#         # proc.plot_m0_MALM(k_MALM[i],nodes,imin,m0[j])


# %%

# sol = []
# for i, f in enumerate(filenames_MALM):
#     sol.append(proc.invert_MALM(imaging_path + 'inversionMALM/' + f,
#                                 pareto=True)
#                 )
#     proc.plot_MALM(k_MALM[i],nodes,imin,sol[i],show=True)


# %% In[ ]:

# k_TL = proc.invert_ERT_TL(
#                             imaging_path,
#                             files=filenames_ERT,
#                             regType=1,
#                             recip=10
#                         )

# proc.plot_ERT(k_TL, vmin=0, vmax=3, attr="Resistivity(log10)", index=0)
# proc.plot_ERT(k_TL, vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=0)
# proc.plot_ERT(k_TL, vmin=0, vmax=50, attr="Resistivity(ohm.m)", index=1)
# proc.plot_ERT(k_TL, vmin=-50, vmax=50, attr="difference(percent)", index=1)
