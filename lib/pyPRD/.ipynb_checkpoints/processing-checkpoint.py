import warnings
warnings.filterwarnings('ignore')
import os
import sys
sys.path.append((os.path.relpath('../src'))) # add here the relative path of the API folder
testdir = './rawData/'
from resipy import Project
import numpy as np




def invert_PRD(path='/home/ben/Documents/Ongoing/Rhizotron_PRD/PRD/rawData/',
               timeStep=2):

    name = 'PRD' + str(timeStep)
    
    if not os.path.exists('inversion/' + name):
        os.makedirs('inversion/' + name)


    k = Project('inversion/' + name, typ="R3t")
    # k.createSurvey('/home/ben/Documents/Ongoing/Rhizotron_PRD/PRD/rawData/'+ name + '.csv', ftype="Syscal")
    k.createSurvey(path + 'rawData/' + name + '.csv', ftype="Syscal")
    k.showError(index=0)
    k.filterRecip(percent=10)
    # k.showErrorDist()
    # k.filterManual(index=0, attr="resist" # interactive plot)
    k.showPseudo()
    k.importElec(path + "mesh/elecsXYZ.csv")
    #k.setElec(elec)
    # k.filterElec([28])
    k.showPseudo()
    # k.filterManual(index=0, attr="resist" # interactive plot)
    k.importMesh(path + "mesh/mesh_rhizo_resipy.msh")
    # k.showMesh()
    k.setStartingRes(regionValues={1: 100.0}, 
                     zoneValues={1: 1}, 
                     fixedValues={1: False}, 
                     ipValues={1: 0.0}) # define region manually using k.addRegion()
    k.param['num_xy_poly'] = 0 # tells R3t to no crop the mesh
    k.param['zmin'] = -np.inf
    k.param['zmax'] = np.inf

    # k.err = True # if we want to use the error from the error models fitted before
    k.param['a_wgt'] = 0.01
    k.param['b_wgt'] = 0.1

    k.invert(modErr=False, parallel=True, modelDOI=False)

    # import matplotlib.pylab as plt
    # fig, ax = plt.subplots()
    import pyvista
    plotter = pyvista.Plotter()
    k.showResults(index=0, attr="Resistivity(log10)",
                  edge_color="none", 
                  vmin=None, 
                  vmax=None, 
                  color_map="viridis", 
                  background_color=(0.8, 0.8, 0.8),
                  pvslices=([], [], []), 
                  pvthreshold=None, 
                  pvdelaunay3d=False, 
                  pvgrid=False, 
                  pvcontour=[],
                  ax=plotter)
    # plotter.save_graphic("img.svg")  

    k.showInvError(index=0)

    # k.saveMeshVtk()
