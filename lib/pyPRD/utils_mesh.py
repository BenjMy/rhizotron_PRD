#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:14:52 2022

@author: ben
"""

# Point(	72	)	=	{	0.4	,	0	,	0.075	,	cl1	};	

imin, nodes, grid = proc.create_source_grid(mesh, k_MALM,
                                            nVRTe_rows=9, nVRTe_cols=9,)

text_file = open("sample.txt", "wt")

for i, g in enumerate(grid):
    pt_nb = i + 72
    defaut_line = 'Point( {} ) = {{{}, {}, {}, cl1}}; \n'.format(pt_nb, g[0],g[1],g[2])
    n = text_file.write(defaut_line)


    
text_file.close()

stri = ''
for i in range(90+72):
    print(i)
    if i==(90+72):
        stri += str(i)
    else:
        stri += str(i) + ','
