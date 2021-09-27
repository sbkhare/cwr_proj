# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:00:53 2020

@author: Sikander
"""

import cwr_model as cwr
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pickle as pkl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl


SINGLE_RUN = 1
BETA_DISTR = 0
SURVIVALPROB = 0
P_SURFACE_PLOT = 0


def F(n, Q):
    return 1 - Q**n

def calc_survival_prob(b0_res, b0_mut, d, mu, save=False):
    N_list = np.arange(500, 10500, 2000, dtype=int)
    K_list = b0_mut*N_list/(b0_mut - d)
    num_trials = 1000
    survival_prob = {} #pkl.dump
    for K, N in zip(K_list, N_list):
        print(">>>    K = ", K, ", N* = ", N)
        Kvec = K*np.ones(num_trials, dtype=float) 
        comm = cwr.community(Kvec, timesteps=200, change="b", dd="b", r=d, res=b0_res, mut=b0_mut, mu=mu)
        comm.simulate_ddb()
        nonzero = np.count_nonzero(comm.N_t[:,-1])
        survival_prob[np.round((1-d/b0_mut)*K)] = nonzero/num_trials #put in N instead of np.round((1-d/b0_mut)*K)?
    
    plt.figure()
    plt.scatter(list(survival_prob.keys()), list(survival_prob.values()), color='k')
    plt.xlabel("Starting population, n")
    plt.ylabel("Rescue frequency")
    
    popt, pcov = curve_fit(F, list(survival_prob.keys()), list(survival_prob.values()))
    n_list = (1-d/b0_mut)*K_list #np.arange(10, 750, 10, dtype=int)
    plt.plot(n_list, F(n_list, *popt), 'r-')
    P = 1-popt[0]
    plt.title("$b_{0,m}$" + " = {0}, d = {1}, P = {2}".format(round(b0_mut, 1), round(d, 1), round(P, 5))) #Rescue probability = 1 - $(1 - P)^n$, where
    plt.grid()
    if save:
        plt.savefig("survival_prob/bm{0}_d{1}.png".format(round(b0_mut, 1), round(d, 1)))
        plt.close()
    return P

if SINGLE_RUN == 1:
    b0_r = 0.05
    b0_m = 0.6 #0.2, 0.4, 0.6, 0.8
    d_all = 1-np.exp(-0.1)
    mu_all = 1-np.exp(-0.0002)
    
    num_species = 167
    
    K_global = 16936#750 # #106 #
    
    #Same K
#    N_vec = K_global*np.ones(num_species, dtype=int) 
    
#    Lognormal
    N_vec = np.random.lognormal(3, 2, num_species) #5, 2
    
    #Beta
    # N_vec = K_global*np.random.beta(0.2498, 41.4668, num_species)
    
    #Turn N_vec into K_vec
    K_vec = b0_m*N_vec/(b0_m - d_all)
    K_vec = K_vec.astype(int)
    
    c = cwr.community(K=K_vec, timesteps=200, change="b", dd="b", r=d_all, res=b0_r, mut=b0_m, mu=mu_all)
    c.simulate_ddb()
    
#     c.plot(save=True)
#     c.plot_diversity(save=True)
#     c.mutation_extinction(save=True)
# #    c.compare_rankabundance()
#     c.species_abundance(0, save=True)
#     c.species_abundance(50, save=True)
#     c.species_abundance(100, save=True)
#     c.species_abundance(150, save=True)
#     c.species_abundance(200, save=True)
# #    c.species_abundance(250)
    c.plot_fig2([0,50,100,150,200], save=True)
    
    
    
if BETA_DISTR == 1: #Only first one applies to K_global ~ 16000
    K_global = 16000
    
    alphas = [0.2498,0.3868,0.2783,0.4872,0.4291,0.2773,0.2538,0.3025,0.3709,
              0.8796,0.9877,0.6976,0.6427,0.8557,1.1268]
    betas = [41.4668,261.8370,85.4514,399.4604,430.3599,62.1201,62.4323,74.4140,
             41.9143,214.6725,493.8308,124.8712,288.5661,344.0007,300.8603]
    
#    for alpha, beta in zip(alphas, betas):
    
if SURVIVALPROB == 1:
    b0_r = 0.05
    b0_m = 0.6 #0.2, 0.4, 0.6, 0.8
    d_all = 0.1
    mu_all = 0.0002
    
#    K_list = np.arange(25, 800, 25, dtype=int)
#    num_trials = 1000
#    survival_prob = {} #pkl.dump
#    for K in K_list:
#        print(">>>    K = ", K)
#        Kvec = K*np.ones(num_trials, dtype=int) 
#        comm = cwr.community(Kvec, timesteps=100, change="b", dd="b", r=d, res=b0_res, mut=b0_mut, mu=mu)
#        comm.simulate_ddb()
#        nonzero = np.count_nonzero(comm.N_t[:,-1])
#        survival_prob[np.round((1-d/b0_mut)*K)] = nonzero/num_trials
#    
#    plt.figure()
#    plt.scatter(list(survival_prob.keys()), list(survival_prob.values()), color='k')
#    plt.xlabel("Starting population, n")
#    plt.ylabel("Rescue frequency")
#    
#    popt, pcov = curve_fit(F, list(survival_prob.keys()), list(survival_prob.values()))
#    n_list = (1-d/b0_mut)*K_list #np.arange(10, 750, 10, dtype=int)
#    plt.plot(n_list, F(n_list, *popt), 'r-')
#    P = 1-popt[0]
#    plt.title("Rescue probability = 1 - $(1 - P)^n$, where P = " + str(round(P, 4)))
#    plt.grid()
    
    calc_survival_prob(b0_r, b0_m, d_all, mu_all, save=True)
    
if P_SURFACE_PLOT == 1:
    b0_r = 0.05
    mu_ = 0.0002
    bm_list = np.arange(0.06, 1.01, 0.05)
    d_list = np.arange(0.05, 1, 0.05)
    D, B = np.meshgrid(d_list, bm_list)
    # P_grid = np.zeros((len(d_list), len(bm_list)), dtype=float)
    
    
    # for i, bm in enumerate(bm_list):
    #     for j, d_ in enumerate(d_list):
    #         print("b0_mut = {0}, d = {1}".format(bm, d_))
    #         if d_ >= bm:
    #             P_grid[i,j] = 0
    #         else:
    #             P = calc_survival_prob(b0_r, bm, d_, mu_, save=True)
    #             P_grid[i,j] = P
    
    
    fig = plt.figure(figsize=(7, 3), dpi=80)
    u_list = np.arange(0.0001, 0.0009, 0.0002)
    br_list = np.arange(0.01, 0.05, 0.01)
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
        ax1.plot_surface(D, B, P_grid, label="${\mu}$ = " + str(u)) # c = 
        # c._facecolors2d=c._facecolors3d
        # c._edgecolors2d=c._edgecolors3d
    ax1.set_xlabel('d', fontsize=16)
    ax1.set_ylabel("$b_{0,m}$", fontsize=16)
    ax1.set_zlabel('P', fontsize=16)
    # ax1.legend()
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
        ax2.plot_surface(D, B, P_grid, label="$b_{r}$ = " + str(br)) #c = 
        # c._facecolors2d=c._facecolors3d
        # c._edgecolors2d=c._edgecolors3d
    ax2.set_xlabel('d', fontsize=16)
    ax2.set_ylabel("$b_{0,m}$", fontsize=16)
    ax2.set_zlabel('P', fontsize=16)    
    # ax2.legend()
    plt.tight_layout()
    # plt.show()
    # plt.savefig("Psurface.png")
    
    
    
    
    
    
    
    
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    