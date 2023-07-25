#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:17:17 2023

@author: ben
"""

from resipy import Project
# import pygimli as pg

MESH = "../imaging_ERT_MALM/mesh/mesh_rhizo_resipy_vrte_v7.msh"

import resipy.meshTools as mt

mesh = mt.readMesh(MESH)

mesh['node_x']
c = (mesh.elmCentre)
len(c)

41878
len(mesh.cells)

mesh = Mesh(node_x = mesh_dict['node_x'],
            node_y = mesh_dict['node_y'],
            node_z = mesh_dict['node_z'],
            node_data = np.array(mesh_dict['node_data']).T,
            cell_type = mesh_dict['cell_type']
            