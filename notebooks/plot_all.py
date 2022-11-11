#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 13:27:29 2022

@author: ben
"""


import imageio
from pathlib import Path
import matplotlib.pyplot as plt
from natsort import os_sorted
import surveyPRD
import argparse
import numpy as np
from PIL import Image
import pickle


def get_cmd():
    parse = argparse.ArgumentParser()
    
    period = parse.add_argument_group('period')
    # period.add_argument('-cycle', type=list, help='cycle to analyse', default=[None], required=False)
    period.add_argument('-cycle','--cycle', nargs='+', help='list of cycle', type=int, default='5', required=False)
    period.add_argument('-startD', type=str, help='start date to analyse', default=None, required=False)
    period.add_argument('-endD', type=str, help='end date', default=None, required=False)
    args = parse.parse_args()

    return(args)

args = get_cmd()

f_ERT, f_MALM,  ERT_log, irr_log = surveyPRD.load_log_files(args)
#%%


#%%
im_fig = []
im_ax = []
pathfile = []
for f in f_MALM:
    # for path in Path('../imaging_ERT_MALM/inversionMALM').rglob(f+ '.png'):
    for path in Path('../imaging_ERT_MALM/inversionMALM').rglob(f+ '.obj'):
        pathfile.append(path)
    
for p in os_sorted(pathfile):
    # im.append(imageio.imread(p))
    with open(p, 'rb') as file:
        im_fig.append(pickle.load(file)[0])
    with open(p, 'rb') as file:
        im_ax.append(pickle.load(file)[0])
    
#%%

# fig = plt.figure(constrained_layout=True)
# ax_array = fig.subplots(2, 2, squeeze=False)

# ax_array[0, 0].bar(["a", "b", "c"], [5, 7, 9])
# ax_array[0, 1].plot([1, 2, 3])
# ax_array[1, 1].imshow([[1, 2], [2, 1]])





letters = 'abcdefghijkl'
l_ax = []
for i in range(len(im_ax)):
    l_ax.append([letters[i]])

fig_handle = pickle.load(open(pathfile[0], "rb"))

fig = plt.figure()
ax_master = fig_handle[0]

type(fig)

fig_handle[0].show()
# axs = plt.subplots([ax_master,ax_master])

fig_new, axs = plt.subplot_mosaic(l_ax,
                              constrained_layout=True)
plt.sca(fig_handle[0].axes[0])        
        # set the current axes instance to the top left


# axs[1][1] = ax_master 
# for ax in fig_handle[0].axes:
#     if ax is not ax_master:
#         ax_master.get_shared_y_axes().join(ax_master, ax)

fig.show()



# ax1[0].gca()

# fig.canvas.draw()

fig_handle[0].axes[0].lines[0].get_data()
fig_handle[0].axes[0].images[0].get_data()

# fig, axs = plt.subplot_mosaic(l_ax,
#                               constrained_layout=True)
# type(ax1[0])
for i, imid in enumerate(im_ax):
    # axs[l_ax[i][0]]= ax1[1]
    axs[1][i]= ax1[1]
    # axs[1][i].plot(ax1[0])
fig.canvas.draw()

# ax2 = plt.axes()




# plt.show()
# plt.savefig('../imaging_ERT_MALM/inversionMALM/myimage.png', dpi=450)


#%%

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
np.arange(4,len(im))


fig, axs = plt.subplot_mosaic([['a)', 'c)'], ['b)', 'c)'], ['d)', 'd)']],
                              constrained_layout=True)

for label, ax in axs.items():
    # label physical distance in and down:
    trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
            fontsize='medium', verticalalignment='top', fontfamily='serif',
            bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

plt.show()



#%%
im = []
pathfile = []
for path in Path('../imaging_ERT_MALM/inversionERT').rglob("Resistivity(ohm.m)_t0.png"):
    pathfile.append(path)
    
for p in os_sorted(pathfile):
    im.append(imageio.imread(p))


plt.figure()
plt.subplots_adjust(right=-0.2)

for i, j in enumerate(np.arange(4,len(im))):
    ax = plt.subplot(3, 2, i + 1)
    # plt.subplot(5,5,i+1)    # the number of images in the grid is 5*5 (25)
    plt.imshow(im[i])
    plt.axis('off')
    
plt.show()
# plt.savefig('test.png')
plt.savefig("myimage.eps", dpi=1200)
plt.savefig("myimage.png", dpi=1200)

