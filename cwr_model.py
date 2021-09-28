# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 19:41:41 2020

@author: Sikander
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

 # 0.00005, 0.005, 0.5
        
class community():
    def __init__(self, K, timesteps, change, dd, r, res, mut, mu):
        self.K = K #vector of species K
        self.timesteps = timesteps
        self.species_vector = np.zeros(len(K), dtype=object) #array of lists
        self.N_t = np.zeros((len(K), timesteps+1), dtype=int) #records species populations over time
        self.H_t = np.zeros(timesteps+1, dtype=float) #records species diversity over time
        self.J_t = np.zeros(timesteps+1, dtype=float) #records species evennes accross time
        self.first_mut = {} #records when first mutation arises
        self.extinction_time = {} #records time to extinction SO THAT MEDIAN TIME CAN BE CALCULATED
        self.dd = dd
        self.mu = mu
        self.change = change #does environmental change affect b or d
        if change == "b":
            self.d = r
            self.b0_res = res
            self.b0_mut = mut
            self.N_t[:,0] = np.round((1-self.d/self.b0_mut)*self.K)
            H, J = self.diversity_evenness(0)
            self.H_t[0] = H
            self.J_t[0] = J
        elif change == "d":
            self.b = r
            self.d0_res = res
            self.d0_mut = mut
        
            
        
        
    def step_ddb(self, stepnumber): #ddb = density dependence on birth
        # print("t = " + str(stepnumber))
        # new_species_vector = np.zeros(len(self.K), dtype=object)
        # for i, species in enumerate(self.species_vector):
        #     N_i = len(species)
        #     new_indiv = []
        #     dead_indiv = []
        #     if i not in self.extinction_time:
        #         if len(species) == 0:
        #             self.extinction_time[i] = stepnumber - 1
        #     for j, indiv in enumerate(species): #iterates through the specieis individual list
        #         b_i = max(indiv*(1 - N_i/self.K[i]), 0) #birth rate can't be lower than 0
        #         if np.random.rand() < b_i: #if random number is less than b0
        #             if indiv == self.b0_res and np.random.rand() < self.mu: #mutation conditional on birth of a resident
        #                 new_indiv.append(self.b0_mut)
        #                 if i not in self.first_mut: #If species index isnt in the dct, add it with stepnumvber as the value
        #                     self.first_mut[i] = stepnumber
        #             else:
        #                 new_indiv.append(indiv)
        #         if np.random.rand() < self.d:
        #             dead_indiv.append(j) #records list index of dead individual
        #     new_species_vector[i] = self.species_vector[i] + new_indiv
        #     for index in sorted(dead_indiv, reverse=True): #Pops dead individuals from list in descending order of index so that lower indices are not changed after pop
        #         new_species_vector[i].pop(index)
        # self.species_vector = new_species_vector
        # for i, species in enumerate(self.species_vector):
        #     self.N_t[i, stepnumber] = len(species) #Records species poulations at stepnumber
        # H, J = self.diversity_evenness(stepnumber)
        # self.H_t[stepnumber] = H
        # self.J_t[stepnumber] = J
        
        print("t = " + str(stepnumber))
        new_species_vector = np.zeros(len(self.K), dtype=object)
        for i, species in enumerate(self.species_vector):
            N_i = len(species)
            new_indiv = []
            dead_indiv = []
            mutated_indiv = []
            if i not in self.extinction_time:
                if len(species) == 0:
                    self.extinction_time[i] = stepnumber - 1
            for j, indiv in enumerate(species): #iterates through the specieis individual list
                b_i = indiv*(1 - N_i/self.K[i])
                if np.random.rand() < self.mu:
                    mutated_indiv.append(j)
                    if i not in self.first_mut: #If species index isnt in the dct, add it with stepnumvber as the value
                        self.first_mut[i] = stepnumber
                if np.random.rand() < b_i: #if random number is less than b0
                    new_indiv.append(indiv)
                if np.random.rand() < self.d:
                    dead_indiv.append(j) #records list index of dead individual
            new_species_vector[i] = self.species_vector[i] + new_indiv
            for index in mutated_indiv: #Replaces individual if muated
                new_species_vector[i][index] = self.b0_mut
            for index in sorted(dead_indiv, reverse=True): #Pops dead individuals from list in descending order of index so that lower indices are not changed after pop
                new_species_vector[i].pop(index)
        self.species_vector = new_species_vector
        for i, species in enumerate(self.species_vector):
            self.N_t[i, stepnumber] = len(species) #Records species poulations at stepnumber
        H, J = self.diversity_evenness(stepnumber)
        self.H_t[stepnumber] = H
        self.J_t[stepnumber] = J
    
    def simulate_ddb(self):
        print("Total community size at t=0: ", np.sum(self.N_t[:,0]))
        for i, species in enumerate(list(self.species_vector)):
            self.species_vector[i] = self.N_t[i,0]*[self.b0_res]
        for i in range(1, self.timesteps+1):
            self.step_ddb(i)
        print("Total community size at the end: ", np.sum(self.N_t[:,-1]))
            
    def plot(self, save=False):
        plt.figure()
        plt.xlabel("Timestep")
        plt.ylabel("Species abundance")
        for species in self.N_t:
            plt.plot(np.arange(self.timesteps+1), species)
#        for idx in self.first_mut:
#            plt.axvline(x=self.first_mut[idx], color='r', linestyle='dashed', label=idx)
        if save:
            plt.savefig("single_run/All_N_t.png")
    
    def RAC(self, times, panel, save=False):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        rank = np.arange(1, len(self.K)+1)
        ax.text(-0.1, 1.15, panel, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        for t in times:
            abundance = self.N_t[:,t]
            ax.plot(rank, sorted(abundance, reverse=True), linewidth=5, label="t = {0}".format(t))
            ax.set_xlabel("Rank")
            ax.set_ylabel("Abundance")
            ax.set_yscale("log")
        ax.legend()
        #Fit beta distribution here
        if save:
            plt.savefig("single_run/rac.png")
        
    def species_abundance(self, t, save=False):
        ###NEEED TO REWRITE############
        plt.figure()
        N_0 = self.N_t[:,0]
        rank_0 = sorted(range(len(N_0)), reverse=True, key=lambda k: N_0[k]) #sorts species indices by descending order of abundance at t=0
        N_t = []
        for index in rank_0:
            N_t.append(self.N_t[index,t])
        plt.bar(range(len(rank_0)), N_t) #For some reason, it automatically sorts x labels even when they are strings, so I have to add xtixks manually
#        plt.xticks(range(len(rank_0)), rank_0)
        plt.title("Species abundance, t=" + str(t))
        plt.xlabel("Species indexed according to initial rank")
        plt.ylabel("Abundance")
        plt.yscale("log")
        if save:
            plt.savefig("single_run/sac{0}.png".format(t))
        
    def diversity_evenness(self, t):
        """https://en.wikipedia.org/wiki/Diversity_index#Shannon_index
        https://en.wikipedia.org/wiki/Species_evenness"""
        N_t = []
        for Ni in self.N_t[:,t]:
            if Ni > 0:
                N_t.append(Ni)
        N_t = np.array(N_t)
        total = np.sum(N_t)
        p_t = N_t/total #NEED TO REMOVE ZEROS
        H = -1*np.sum(p_t*np.log(p_t)) #diversity
        S = len(N_t) #number of extant species
        Hmax = np.log(S) 
        J = H/Hmax #eveness 
        return H, J
    
    def plot_diversity(self, save=False):
        plt.figure()
        plt.subplot(121)
        plt.plot(np.arange(self.timesteps+1), self.H_t, 'k')
        plt.xlabel("Timestep")
        plt.ylabel("Shannon's diversity index")
        plt.subplot(122)
        plt.plot(np.arange(self.timesteps+1), self.J_t, 'k')
        plt.xlabel("Timestep")
        plt.ylabel("Pielou's evenness index")
        
        #Calculate rarity gain
        # N_0 = []
        # N_f = []
        # for Ni1, Ni2 in zip(self.N_t[:,0], self.N_t[:,-1]):
        #     if Ni1 > 0 and Ni2 > 0:
        #         N_0.append(Ni1)
        #         N_f.append(Ni2)
        # N_0 = np.array(N_0)
        # rank_0 = sorted(range(len(N_0)), reverse=True, key=lambda k: N_0[k]) #sorts species indices by descending order of abundance at t=0
        # N_f = np.array(N_f)
        # tot_0 = np.sum(N_0)
        # tot_f = np.sum(N_f)
        # p_0 = N_0/tot_0
        # p_f = N_f/tot_f
        # rarity_change = np.log(p_f) - np.log(p_0)
        # Tau = np.sum(p_f*rarity_change)
        # del_rar = []
        # for index in rank_0:
        #     del_rar.append(rarity_change[index])
        # plt.subplot(133)
        # plt.bar(range(len(rank_0)), del_rar)
        # plt.title("Rarity/information gain = " + str(Tau))
        # plt.xlabel("Species ranked according to initial abundance")
        # plt.ylabel("Change in rarity")
        
        if save:
            plt.savefig("single_run/diversity_evenness.png")
        
    def mutation_extinction(self, save=False):
        print("Median time to first mutation = ", np.median(list(self.first_mut.values())))
        print("Median extinction time = ", np.median(list(self.extinction_time.values())))
        initial_pop = []
        for idx in self.first_mut:
            initial_pop.append(self.N_t[idx,0])
        plt.figure()
        plt.subplot(121)
        plt.scatter(initial_pop, self.first_mut.values(), color='blue')
        plt.xlabel("Initial population")
        plt.ylabel("Time of first mutation")
        plt.subplot(122)
        initial_pop = []
        for idx in self.extinction_time:
            initial_pop.append(self.N_t[idx,0])
        plt.scatter(initial_pop, self.extinction_time.values(), color='green')
        plt.xlabel("Initial population")
        plt.ylabel("Time to extinction")
        if save:
            plt.savefig("single_run/mutation_extinction.png")
        plt.figure()
        ext_times = list(self.extinction_time.values())
        ext_times = list(filter(lambda a: a != 0, ext_times))
        plt.hist(ext_times, color='orange', bins=np.arange(0, self.timesteps, 5))
        plt.xlabel("Time")
        plt.ylabel("Number of extinctions")
        if save:
            plt.savefig("single_run/extinction_debt.png")
            
    def extinction_debt(self, title, panel, col, save=False):
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.text(-0.1, 1.15, panel, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        # ax.set_title(title)
        ext_times = list(self.extinction_time.values())
        ext_times = list(filter(lambda a: a != 0, ext_times))
        ax.hist(ext_times, color=col, bins=np.arange(0, self.timesteps, 5))
        plt.xlabel("Time")
        plt.ylabel("Number of extinctions")
        if save:
            plt.savefig("single_run/extinction_debt_{0}.png".format(title))
        
#    def compare_rankabundance(self):
#        N_0 = []
#        for Ni in self.N_t[:,0]:
#            if Ni > 0:
#                N_0.append(Ni)
#        rank_0 = sorted(range(len(N_0)), reverse=True, key=lambda k: N_0[k])
#        N_end = []
#        for Ni in self.N_t[:,-1]:
#            if Ni > 0:
#                N_end.append(Ni)
#        rank_end = sorted(range(len(N_end)), reverse=True, key=lambda k: N_end[k])
#        rank_dict = {} #Keys:initial rank, vals: end rank
#        for rend, idx in enumerate(rank_end):
#            r0 = rank_0.index(idx)
#            rank_dict[r0] = rend
#        plt.figure()
#        plt.scatter(rank_dict.keys(), rank_dict.values())
#        plt.xlabel("Initial rank")
#        plt.ylabel("Final rank")
#        res = stats.linregress(list(rank_dict.keys()), list(rank_dict.values()))
#        x = np.linspace(min(rank_dict.keys()), max(rank_dict.keys()), 1)
#        plt.plot(x, res.intercept + res.slope*x, 'r', label='fitted line')
#        plt.legend()
#        plt.title("R-squared: {0}, p-value: {1}".format(round(res.rvalue**2, 3), res.pvalue**2))

    def plot_fig2(self, times, save=False): 
        fig = plt.figure(figsize=(7, 6), dpi=80)
        ax1 = fig.add_subplot(2, 2, 1)
        ax1.text(-0.1, 1.15, "a", transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        ax1.set_xlabel("Timestep")
        ax1.set_ylabel("Species abundance")
        for species in self.N_t:
            ax1.plot(np.arange(self.timesteps+1), species)
        ax2 = fig.add_subplot(2, 2, 2)
        ax2.text(-0.1, 1.15, "b", transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        rank = np.arange(1, len(self.K)+1)
        for t in times:
            abundance = self.N_t[:,t]
            ax2.plot(rank, sorted(abundance, reverse=True), linewidth=5, label="t = {0}".format(t))
            ax2.set_xlabel("Rank")
            ax2.set_ylabel("Abundance")
            ax2.set_yscale("log")
        ax2.legend()
        ax3 = fig.add_subplot(2, 2, 3)
        ax3.text(-0.1, 1.15, "c", transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        ax3.plot(np.arange(self.timesteps+1), self.H_t, 'k')
        ax3.set_xlabel("Timestep")
        ax3.set_ylabel("Diversity")
        ax4 = fig.add_subplot(2, 2, 4)
        ax4.text(-0.1, 1.15, "d", transform=ax4.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        ax4.plot(np.arange(self.timesteps+1), self.J_t, 'k')
        ax4.set_xlabel("Timestep")
        ax4.set_ylabel("Evenness")
        plt.tight_layout()
        if save:
            plt.savefig("single_run/fig2.png")
        
                
if __name__=='__main__':
    
    ##########      Run single simulation       ###############################
    num_species = 150
    K_global = 16000
    
    
#    K_vec = K_global*np.ones(num_species, dtype=int) 
#    
    K_vec = np.random.lognormal(3, 1, num_species)
    K_vec = K_global*np.random.beta(0.2498, 41.4668)
    K_vec = K_vec.astype(int)
#    
#    K_vec = 
    
    c = community(K_vec, timesteps=100)
    c.simulate()
    c.plot()
#    c.RAC(0)
#    c.RAC(25)
    c.species_abundance(0)
    c.species_abundance(25)
    c.species_abundance(50)
    c.species_abundance(75)
    c.species_abundance(100)
    
    
    
    ########    Plot K vs survival probaiblity    #############################
    K_list = np.arange(10, 750, 10, dtype=int)
    num_trials = 1000
    survival_prob = {}
    for K in K_list:
        K_vec = K*np.ones(num_trials, dtype=int) 
        c = community(K_vec, timesteps=100)
        c.simulate()
        nonzero = np.count_nonzero(c.N_t[:,-1])
        survival_prob[K] = nonzero/num_trials
    plt.figure()
    plt.scatter(list(survival_prob.keys()), list(survival_prob.values()), color='k')
#    xp = np.linspace(5,255, 5)
#    m, b = np.polyfit(list(survival_prob.keys()), list(survival_prob.values()), 1)
#    plt.plot(xp, m*xp+b, color='r', linestyle='-')
    plt.xlabel("K")
    plt.ylabel("Probability of rescue")
        
        
    
    
