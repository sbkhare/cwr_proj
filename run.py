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

SINGLE_RUN = 1
BETA_DISTR = 0
SURVIVALPROB = 0
P_SURFACE_PLOT = 0


def F(n, Q):
    return 1 - Q**n

def calc_survival_prob(b0_res, b0_mut, d, mu, save=False):
    N_list = np.arange(1000, 30000, 1000, dtype=int)
    K_list = b0_mut*N_list/(b0_mut - d)
    num_trials = 1000
    survival_prob = {} #pkl.dump
    for K, N in zip(K_list, N_list):
        print(">>>    K = ", K, ", N = ", N)
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
    b0_m = 0.9 #0.2, 0.4, 0.6, 0.8
    d_all = 0.1
    mu_all = 0.0002
    
    num_species = 247
    
    K_global = 23756 #750 # #106 #
    
    #Same K
#    N_vec = K_global*np.ones(num_species, dtype=int) 
    
#    Lognormal
#    N_vec = np.random.lognormal(4.15, 2, num_species) #5, 2
    
    #Beta
    N_vec = K_global*np.random.beta(0.2538, 0.2538, num_species)
    
    #Turn N_vec into K_vec
    K_vec = b0_m*N_vec/(b0_m - d_all)
    K_vec = K_vec.astype(int)
    
    c = cwr.community(K=K_vec, timesteps=200, change="b", dd="b", r=d_all, res=b0_r, mut=b0_m, mu=mu_all)
    c.simulate_ddb()
    c.plot(save=True)
    c.plot_diversity(save=True)
    c.mutation_extinction(save=True)
#    c.compare_rankabundance()
    c.RAC(0, save=True)
    c.RAC(50, save=True)
    c.RAC(100, save=True)
    c.RAC(150, save=True)
    c.RAC(200, save=True)
#    c.RAC(250)
    c.species_abundance(0, save=True)
    c.species_abundance(50, save=True)
    c.species_abundance(100, save=True)
    c.species_abundance(150, save=True)
    c.species_abundance(200, save=True)
#    c.species_abundance(250)
    
    
    
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
    bm_list = np.arange(0.1, 1, 0.1)
    d_list = np.arange(0.1, 1, 0.1)
    D, B = np.meshgrid(d_list, bm_list)
    P_grid = np.zeros((len(d_list), len(bm_list)), dtype=float)
    for i, bm in enumerate(bm_list):
        for j, d_ in enumerate(d_list):
            print("b0_mut = {0}, d = {1}".format(bm, d_))
            if d_ >= bm:
                P_grid[i,j] = 0
            else:
                P = calc_survival_prob(b0_r, bm, d_, mu_, save=True)
                P_grid[i,j] = P
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(D, B, P_grid) 
    ax.set_xlabel('d')
    ax.set_ylabel("$b_{0,m}$")
    ax.set_zlabel('P')
    plt.savefig("Psurface_br{0}_mu{0}.png".format(b0_r, mu_))
    
    
    
    
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    