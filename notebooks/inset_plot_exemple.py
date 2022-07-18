#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 21:23:32 2022

@author: ben
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

time = np.arange(0,1e2)
value1 = np.sin(time/2)
value2 = np.cos(time/10)+3
value3 = np.cos(time/6)+3


# T, Cp = np.loadtxt('Ta-Cp.txt', unpack=True)
# T_E, CV_E = np.loadtxt('Ta-CV_Einstein.txt', unpack=True)
# T_D, CV_D = np.loadtxt('Ta-CV_Debye.txt', unpack=True)

fig, ax1 = plt.subplots()
# T_E = np.arange(1,max(T)+1,1)
# The data.
ax1.plot(time, value1, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
# The Einstein fit.
ax1.plot(time, value1, c='m', lw=2, alpha=0.5, label='Einstein model')
ax1.set_xlabel(r'$T\,/\mathrm{K}$')
ax1.set_ylabel(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')
ax1.legend(loc=0)

# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
# ax2 = plt.axes([0,0,1,1])
# # Manually set the position and relative size of the inset axes within ax1
# ip = InsetPosition(ax1, [0.4,0.2,0.5,0.5])
# ax2.set_axes_locator(ip)
# # Mark the region corresponding to the inset axes on ax1 and draw lines
# # in grey linking the two axes.
# mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')


ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(ax1, [0.2,0.4,0.5,0.5])
ax2.set_axes_locator(ip)

mark_inset(ax1, ax2, loc1=2, loc2=4, fc='red', alpha=0.1, ec='0.1')


# The data: only display for low temperature in the inset figure.
Tmax = max(time/3)
ax2.plot(time[time<=Tmax], value1[time<=Tmax], 'x', c='b', mew=2, alpha=0.8,
          label='Experiment')
# The Einstein fit (not very good at low T).
ax2.plot(time[time<=Tmax], value2[time<=Tmax], c='m', lw=2, alpha=0.5,
          label='Einstein model')
# The Debye fit.
ax2.plot(time/2, value3, c='r', lw=2, alpha=0.5, label='Debye model')
ax2.legend(loc=0)

# Some ad hoc tweaks.
ax1.set_ylim(-1,26)
ax2.set_yticks(np.arange(0,2,0.4))
ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
ax2.tick_params(axis='x', which='major', pad=8)

plt.show()




