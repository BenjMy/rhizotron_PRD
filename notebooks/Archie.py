#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 16:21:44 2022

@author: ben
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# Load lab exp. data
path = '../Archie/'


sname = [str(s) + 'A' for s in np.arange(7)+1]
saturation= [72.8,72.8,61.1,54.1,51.3,45.8,40.7]
saturation01 = [s/100 for s in saturation]

porosity = 0.55
theta = [s*porosity for s in saturation]

Archie_data = pd.DataFrame(np.array([sname,saturation01,theta]).T,
                           columns=['sname','saturation','theta']
                           )
meanRes = []
for s in sname:
    meas = pd.read_excel(path + 'sandpeatwater_all_meas.xlsx', s,nrows=118)
    frequencies = pd.to_numeric(meas[meas.columns[0]])
    
    minfreq = np.where(frequencies==frequencies.min())[0]
    meanRes.append(meas['resistivity [ohm m]'].iloc[minfreq].mean())
    
Archie_data['mean resistivity [ohm m]'] = meanRes

Archie_data = Archie_data.drop(3)

# 594 microS/cm
sw_conductivity = 594*(1e-6/1e-2)
rho_sat = 1/sw_conductivity
# rho_sat0 = sw_conductivity


from scipy.optimize import least_squares


# ρ=aρflϕ−mS−n

def calc_rho(params, sat, rho, porosity, rho_sat):
    m, n = params
    # rho_calc = rho_sat * (1 / (sat ** n))
    # m = 1.5
    rho_calc = rho_sat * porosity**(-m) * sat **(-n)
    misfit = rho - rho_calc
    return(misfit)


def fit_archie(sat, rho, porosity, rho_sat, m0=None, n0=None):
    if m0 is None:
        m = 1.3
    if n0 is None:
        n = 1
    p0 = (m,n)
    bounds = [(1.3,1), (2.5,3)]
    solution = least_squares(
        fun=calc_rho,
        x0=p0,
        jac='3-point',
        bounds=bounds,
        args=(sat, rho, porosity, rho_sat),
        loss='soft_l1',
        verbose=1,
    )
    return(solution)

Archie_data['saturation'] = Archie_data['saturation'].astype('float')
Archie_data['mean resistivity [ohm m]'] = Archie_data['mean resistivity [ohm m]'].astype('float')
Archie_data = Archie_data.sort_values(by='saturation',ascending=True)

fit_sol = fit_archie(Archie_data['saturation'].to_numpy().astype('float'),
                     Archie_data['mean resistivity [ohm m]'].to_numpy().astype('float'),
                     porosity = 0.55,
                     rho_sat=rho_sat
                     )

# Archie_data['mean_res_fit']= np.flip(fit_sol.fun)
Archie_data['residuals_mean_res_fit']= fit_sol.fun

fig, ax = plt.subplots(1)
# Archie_data.plot(x='saturation',y='mean_res_fit',ax=ax)
Archie_data.plot.scatter(x='saturation',y='mean resistivity [ohm m]',ax=ax)

sat_synth =np.linspace(0.3, 1)
m, n = fit_sol.x
print(m,n)
rho_calc_synth = rho_sat * porosity**(-m) * sat_synth **(-n)


# fig, ax = plt.subplots(1)
ax.plot(sat_synth,rho_calc_synth)


def RMS(y_fit,y):
    RMS = np.sqrt(sum(Archie_data['residuals_mean_res_fit'].abs()**2))
    
    # residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)
    
    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    
    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    
    return RMS, r2



y_fit = rho_sat * porosity**(-m) * Archie_data['saturation']  **(-n)

RMS,r2 = RMS(y_fit,Archie_data['mean resistivity [ohm m]'])


R2_text = '$R^{2}$=' + str(np.round(r2,2))

ax.text(0.8,200,R2_text)
ax.set_ylabel('resistivity [$\Omega$.m]')
ax.set_xlabel('Saturation [-]')

# Archie_data = Archie_data.sort_values(by='saturation',ascending=True)
# Archie_data.plot(x='theta',y='mean_res_fit',ax=ax)
# Archie_data.plot.scatter(x='theta',y='mean resistivity [ohm m]',ax=ax)


