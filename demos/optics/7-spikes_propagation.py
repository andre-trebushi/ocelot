#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 14:15:40 2021

@author: andrei
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sc
from scipy import signal, misc
from ocelot.optics.wave import *
from ocelot.gui.dfl_plot import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and

ocelog.setLevel(logging.INFO)

n_s = 200 # number of periods in an undulator
l_w = 0.018 # [m] undulator period 
L_w = l_w * n_s # undulator length

E_ph = 2167 # resonance energy
w = E_ph / hr_eV_s 
xlamds = 2 * np.pi * speed_of_light / w

sigma_r = np.sqrt(2*xlamds*L_w)/4/np.pi #natural radiation size in the waist
sigma_rp = np.sqrt(xlamds/2/L_w) #natural radiation divergence at the waist

ebeam_sigma_x = 38e-06 # x electron beam size 
ebeam_sigma_y = 4.68e-06 # y electron beam size
ebeam_sigma_xp = 25e-06 # x electron beam divergence
ebeam_sigma_yp = 20e-06 # y electron beam divergence

ebeam_sigma_z = 1 # Electron beam duration. Not physical for the current version of
                  # SERVAL as here is considered the case of fully monochromatic radiation
                  # long beam duration corresponds to an infinitesimally small spectral line

#%%
N_b = 1 #number of statistical realizations
N_e = 150 #number of macro electrons 
Nz, Ny, Nx = N_b, 701, 701 # the shape of the dfl.fld

e_beam_param = r'$N_x$ = {}, '.format(round((ebeam_sigma_x)**2/xlamds/L_w, 3)) + r'$N_y$ = {}, '.format(round((ebeam_sigma_y)**2/xlamds/L_w, 3)) + \
               r'$D_x$ = {}, '.format(round((ebeam_sigma_xp)**2 * L_w/xlamds, 3)) + r'$D_y$ = {}, '.format(round((ebeam_sigma_yp)**2 * L_w/xlamds, 3)) + \
               r'$N_b$ = {} '.format(N_b) + r'$N_e = {}$'.format(N_e)

print(e_beam_param)
#%%
Lz, Ly, Lx = 1000e-6, 1500e-5, 1500e-5 #size of realspace grid [m]
dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz

### creating RadiationField object
dfl_50 = RadiationField((Nz, Ny, Nx))
dfl_50.dx, dfl_50.dy, dfl_50.dz = dx, dy, dz
dfl_50.xlamds = xlamds
dfl_50.filePath = ''
dfl_50.to_domain('sf')

fieldname_50 = '1-far_zone_50_m'
# approximation = "far_field"
approximation = "near_field"

dfl_50 = undulator_field_dfl_MP(dfl_50, z=50, L_w=L_w, E_ph=1240, N_e=N_e, N_b=N_b, 
                                            sig_x=ebeam_sigma_x, sig_y=ebeam_sigma_y, sig_xp=ebeam_sigma_xp, 
                                            sig_yp=ebeam_sigma_yp, approximation=approximation, mode='incoh', seed=1)
dfl_50.to_domain(domains='sf') 
plot_dfl(dfl_50, domains='sf', phase=True, fig_name=fieldname_50)
plot_dfl(dfl_50, domains='kf', phase=True, fig_name =fieldname_50)

#%%
# Lz, Ly, Lx = 1e-6, 300e-5 , 200e-5 / 4 #size of realspace grid [m]
Lz, Ly, Lx = 1, 150e-6/4,  150e-6/4 #size of realspace grid [m]

dx, dy, dz = Lx / Nx, Ly / Ny, Lz / Nz

### creating RadiationField object
dfl_10 = RadiationField((Nz, Ny, Nx))
dfl_10.dx, dfl_10.dy, dfl_10.dz = dx, dy, dz
dfl_10.xlamds = xlamds
dfl_10.filePath = filePath
dfl_10.to_domain('sf')

fieldname_10 = '1-far_zone_10_m'

dfl_10 = undulator_field_dfl_MP(dfl_10, z=10, L_w=L_w, E_ph=1240, N_e=N_e, N_b=N_b, 
                                            sig_x=ebeam_sigma_x, sig_y=ebeam_sigma_y, sig_xp=ebeam_sigma_xp,  
                                            sig_yp=ebeam_sigma_yp, approximation=approximation, mode='incoh', seed=1)

dfl_10.to_domain(domains='sf') 
plot_dfl(dfl_10, domains='sf', phase=True, fig_name = fieldname_10)
plot_dfl(dfl_10, domains='kf', phase=True, fig_name = fieldname_10)

#%%
dfl_10_40 = deepcopy(dfl_10)
dfl_10_40.prop_m(40, m=[3,3])
# dfl_10_40.prop(40)

fieldname = '1-far_zone_10+40_m' #+ simulation_name
dfl_10_40.to_domain(domains='sf') 

plot_dfl(dfl_10_40, domains='sf', phase=True, fig_name = fieldname)
plot_dfl(dfl_10_40, domains='kf', phase=True, fig_name = fieldname)

plot_two_dfls(dfl_10_40, dfl_50, domains='s', fig_name='comparison-far_zone_10+40_m', 
              slice_xy=True, phase=False, label_first='10 + 40', 
              label_second='50', title=e_beam_param, filePath=filePath, savefig=True)





















