# -*- coding: utf-8 -*-
"""
@author: Sikander Khare (s.khare@ufl.edu)
Department of Biology, University of Florida
Gainesville, FL 32608
"""

#PACKAGES
import cwr_model as cwr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#SWITCHES
FIG1_S2_S3 = 1
FIGS1 = 0
FIGS4 = 0



if FIG1_S2_S3 ==1:
    #Produces Figures 1, S2, and S3
    
    
    #Default parameters
    b0_r = 0.05
    b0_m = 0.6 #0.2, 0.4, 0.6, 0.8
    d_all = 0.1 #1-np.exp(-0.1)
    mu_all = 0.0002 #1-np.exp(-0.0002)
    
    num_species = 167
    
    J = 16936 #only for lognormal case
    
    
    #Produce Figure 1
    #Species abundances drawn from a beta distribution
    N_vec1 = J*np.random.beta(0.2498, 41.4668, num_species) #Initial species abundances
    K_vec1 = b0_m*N_vec1/(b0_m - d_all) #K parameters for each species
     
    beta = cwr.community(K=K_vec1, timesteps=200, change="b", dd="b", r=d_all, res=b0_r, mut=b0_m, mu=mu_all) #Intialize community
    beta.simulate() #Simulate
    beta.plot_fig1([0,50,100,150,200])
    
    
    #Produce Figure S2
    #Species abundances drawn from a lognormal distribution
    N_vec2 = np.random.lognormal(1.2706*3, 1.2706, num_species) 
    K_vec2 = b0_m*N_vec2/(b0_m - d_all) 
    
    logn = cwr.community(K=K_vec2, timesteps=200, change="b", dd="b", r=d_all, res=b0_r, mut=b0_m, mu=mu_all) 
    logn.simulate() 
    logn.plot_fig1([0,50,100,150,200])
    
    #Produce Figure S3
    #Temporal expression of the extinction debt
    #Panel a
    beta.extinction_debt("beta", "a", "orange")
    #Panel b
    logn.extinction_debt("lognormal", "b", "orange")
    


if FIGS1 == 1:
    #Produce Figure S1
    #Theoretical probability of evolutoinary rescue
    b0_r = 0.05
    mu_ = 0.0002
    bm_list = np.arange(0.06, 1.01, 0.05)
    d_list = np.arange(0.05, 1, 0.05)
    D, B = np.meshgrid(d_list, bm_list)
    
    fig = plt.figure(figsize=(7, 3), dpi=80)
    u_list = np.arange(0.0001, 0.0009, 0.0002)
    br_list = np.arange(0.01, 0.05, 0.01)
    
    #Panel a
    #Default b_0r and different values of mu
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1.text2D(-0.1, 1.1, s="a", transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    for k, u in enumerate(u_list):
        P_grid = np.zeros((len(d_list), len(bm_list)), dtype=float)
        for i, bm in enumerate(bm_list):
           for j, d_ in enumerate(d_list):
                r0 = b0_r/d_
                r0_star = bm/d_
                p_star = 1 - 1/r0_star
                P = u*r0*p_star/(1 - r0)    
                if P < 0:
                    P = 0
                P_grid[i,j] = P
        ax1.plot_surface(D, B, P_grid)
    ax1.set_xlabel('d', fontsize=16)
    ax1.set_ylabel("$b_{0,m}$", fontsize=16)
    ax1.set_zlabel('P', fontsize=16)
    
    #Panel b
    #Default mu and different values of b_0r
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.text2D(1.5, 1.1, "b", transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    for k, br in enumerate(br_list):
        P_grid = np.zeros((len(d_list), len(bm_list)), dtype=float)
        for i, bm in enumerate(bm_list):
           for j, d_ in enumerate(d_list):
                r0 = br/d_
                r0_star = bm/d_
                p_star = 1 - 1/r0_star
                P = mu_*r0*p_star/(1 - r0)    
                if P < 0:
                    P = 0
                P_grid[i,j] = P
        ax2.plot_surface(D, B, P_grid, label="$b_{r}$ = " + str(br))
    ax2.set_xlabel('d', fontsize=16)
    ax2.set_ylabel("$b_{0,m}$", fontsize=16)
    ax2.set_zlabel('P', fontsize=16)
    
    

if FIGS4 == 1:
    #Simulations with communities larger by an order of magnitude
    #Default parameters
    b0_r = 0.05
    b0_m = 0.6 #0.2, 0.4, 0.6, 0.8
    d_all = 0.1 #1-np.exp(-0.1)
    mu_all = 0.0002 #1-np.exp(-0.0002)
    
    num_species = 167
    
    J = 169360 #order of magnitude larger
    
    #Panel a
    #Beta
    N_vec3 = J*np.random.beta(0.2498, 41.4668, num_species) #Initial species abundances
    K_vec3 = b0_m*N_vec3/(b0_m - d_all) #K parameters for each species
    beta2 = cwr.community(K=K_vec3, timesteps=200, change="b", dd="b", r=d_all, res=b0_r, mut=b0_m, mu=mu_all)
    beta2.simulate()
    beta2.RAC(times=[0,50,100,150,200], panel="a")
    
    #Panel b
    ##Lognormal
    N_vec4 = np.random.lognormal(5, 2, num_species) 
    K_vec4 = b0_m*N_vec4/(b0_m - d_all) 
    
    logn2 = cwr.community(K=K_vec4, timesteps=200, change="b", dd="b", r=d_all, res=b0_r, mut=b0_m, mu=mu_all) 
    logn2.simulate()
    logn2.RAC(times=[0,50,100,150,200], panel="b")