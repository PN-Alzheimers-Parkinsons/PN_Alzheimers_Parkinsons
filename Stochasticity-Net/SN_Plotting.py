# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 10:11:55 2020

@author: brand
"""

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from pylab import plot, ylim, xlim, show, xlabel, ylabel, grid

#from scipy.signal import convolve
#get_ipython().run_line_magic('matplotlib', 'qt')
# Only run this cell once to avoid confusion with directories
# Point this to the directory where HFPN.py is relative to your working directory
cwd = os.getcwd() # Get current working directory
root_folder = os.sep + "PN_Alzheimers_Parkinsons"
# Move to 'utils' from current directory position
sys.path.insert(0, cwd[:(cwd.index(root_folder)+len(root_folder))] + os.sep + "Stochasticity-Net" + os.sep)
from Stochastic_Analysis import Analysis


# In[2]:

analysis = {}

# analysis['6M_SDall10_healthy'] = Analysis.load_from_file('6M_SDall10_healthy')
# analysis['6M_SDall5_healthy'] = Analysis.load_from_file('6M_SDall5_healthy')
# analysis['6M_SDall1_healthy'] = Analysis.load_from_file('6M_SDall1_healthy')
# analysis['6M_SDall15_healthy'] = Analysis.load_from_file('6M_SDall15_healthy')
# analysis['6M_SDall20_healthy'] = Analysis.load_from_file('6M_SDall20_healthy')
# analysis['6M_SDall30_healthy'] = Analysis.load_from_file('6M_SDall30_healthy')
# analysis['6M_SDall40_healthy'] = Analysis.load_from_file('6M_SDall40_healthy')
# analysis['SDall50_healthy2'] = Analysis.load_from_file('SDall50_healthy2')

# analysis['6M_SDall1_aged'] = Analysis.load_from_file('6M_SDall1_aged')
# analysis['6M_SDall5_aged'] = Analysis.load_from_file('6M_SDall5_aged')
# analysis['6M_SDall5_aged2'] = Analysis.load_from_file('6M_SDall5_aged2')
# analysis['6M_SDall10_age'] = Analysis.load_from_file('6M_SDall10_age')
# analysis['6M_SDall10_aged2'] = Analysis.load_from_file('6M_SDall10_aged2')
# analysis['6M_SDall15_aged2'] = Analysis.load_from_file('6M_SDall15_aged2')
# analysis['6M_SDall20_aged'] = Analysis.load_from_file('6M_SDall20_aged')
# analysis['6M_SDall20_aged2'] = Analysis.load_from_file('6M_SDall20_aged2')
# analysis['6M_SDall30_aged'] = Analysis.load_from_file('6M_SDall30_aged')
# analysis['6M_SDall30_aged2'] = Analysis.load_from_file('6M_SDall30_aged2')
# analysis['6M_SDall40_aged'] = Analysis.load_from_file('6M_SDall40_aged')
# analysis['6M_SDall40_aged2'] = Analysis.load_from_file('6M_SDall40_aged2')
# analysis['6M_SDall50_aged'] = Analysis.load_from_file('6M_SDall50_aged')
# analysis['6M_SDall50_aged2'] = Analysis.load_from_file('6M_SDall50_aged2')

#APOE
# analysis['6M_SDall0_apoe'] = Analysis.load_from_file('6M_SDall0_apoe')
# analysis['6M_SDall1_apoe'] = Analysis.load_from_file('6M_SDall1_apoe')
# analysis['6M_SDall1_apoe2'] = Analysis.load_from_file('6M_SDall1_apoe2')
# analysis['6M_SDall5_apoe'] = Analysis.load_from_file('6M_SDall5_apoe')
# analysis['6M_SDall5_apoe2'] = Analysis.load_from_file('6M_SDall5_apoe2')
# analysis['6M_SDall10_apoe'] = Analysis.load_from_file('6M_SDall10_apoe')
# analysis['6M_SDall10_apoe2'] = Analysis.load_from_file('6M_SDall10_apoe2')
# analysis['6M_SDall15_apoe'] = Analysis.load_from_file('6M_SDall15_apoe')
# analysis['6M_SDall15_apoe2'] = Analysis.load_from_file('6M_SDall15_apoe2')
# analysis['6M_SDall20_apoe'] = Analysis.load_from_file('6M_SDall20_apoe')
# analysis['6M_SDall20_apoe2'] = Analysis.load_from_file('6M_SDall20_apoe2')
# analysis['6M_SDall30_apoe'] = Analysis.load_from_file('6M_SDall30_apoe')
# analysis['6M_SDall30_apoe2'] = Analysis.load_from_file('6M_SDall30_apoe2')
# analysis['6M_SDall40_apoe'] = Analysis.load_from_file('6M_SDall40_apoe')
# analysis['6M_SDall40_apoe2'] = Analysis.load_from_file('6M_SDall40_apoe2')
# analysis['6M_SDall50_apoe'] = Analysis.load_from_file('6M_SDall50_apoe')
# analysis['6M_SDall50_apoe2'] = Analysis.load_from_file('6M_SDall50_apoe2')

#Whole_Module_10e6
#Whole_Module_10e6_SD_10_percent
# analysis['LRRK2_mut_100000'] = Analysis.load_from_file('LRRK2_mut_100000')
# analysis['DJ1'] = Analysis.load_from_file('pd_DJ1_ds_smallertimestep')
# analysis['lrrk2'] = Analysis.load_from_file('pd_lrrk2_ds_smallertimestep')
# analysis['vps35'] = Analysis.load_from_file('pd_vps35_ds_smallertimestep')
# analysis['gba1_lrrk2'] = Analysis.load_from_file('pd_gba1_lrrk2_ds_smallertimestep')
# analysis['vps35_lrrk2'] = Analysis.load_from_file('pd_vps35_lrrk2_ds_smallertimestep')
# analysis['gba1_vps35_lrrk2'] = Analysis.load_from_file('pd_gba1_vps35_lrrk2_ds_smallertimestep')
# analysis['all_mutations'] = Analysis.load_from_file('pd_all_mutations_ds_smallertimestep')
# analysis['ApoEchol'] = Analysis.load_from_file('pd_2xApoE_ds_smallertimestep')
# analysis['27OHchol'] = Analysis.load_from_file('pd_2x27OH_ds_smallertimestep')
# analysis['ApoE_lrrk2_gba1'] = Analysis.load_from_file('pd_2xApoE_lrrk2_gba1_ds_smallertimestep')
# analysis['27OH_lrrk2_gba1'] = Analysis.load_from_file('pd_2x27OH_lrrk2_gba1_ds_smallertimestep')
# analysis['NPT'] = Analysis.load_from_file('pd_NPT_ds_smallertimestep')
# analysis['DNL'] = Analysis.load_from_file('pd_DNL_ds_smallertimestep')
# analysis['LAMP2A'] = Analysis.load_from_file('pd_LAMP2A_ds_smallertimestep')

# analysis['Ab6mil'] = Analysis.load_from_file('Ab6mil')
# analysis['Ab6milagedSD20simp'] = Analysis.load_from_file('Ab6milagedSD20simp')
# analysis['test'] = Analysis.load_from_file('test')
# analysis['6M_SD10_aged_AB'] = Analysis.load_from_file('6M_SD10_aged_AB')
# analysis['agedSD10AB'] = Analysis.load_from_file('agedSD10AB')
analysis['CD33Aged6t'] = Analysis.load_from_file('CD33Aged6t')
analysis['healthy6t'] = Analysis.load_from_file('healthy6t')


desired_plotting_steps = 600000
# In[3]:

#brandonadded
def plot_time_evolution(self, place_ids = []):#so the idea is you have the place IDs as a list.
        """Plots time-evolution using mean number of tokens for each place in the net.
        TODO: Only plot evolution of places specified in place_ids if place_ids is not empty.
        """
        #Info On the New Plotting Code Below:
        #places is a dictionary. So you can get the actual places using places.keys()
        #In the zip method below, self.petri_net_model.places.values() was changed to self.petri_net_model.places.keys() because this was more readable.
        
        dict_of_tokens = {} #brandonadded. creates dictionary
        for tokens, place in zip(self.timeseries_mean.T, self.petri_net_model.places.keys()):
            dict_of_tokens["{}".format(place)]=tokens#brandonadded. this line assigns each places list of tokens over all timesteps, to the dictionary dict_of_tokens and assigns the dictionary key using whats inside the square brackets.
    
        for x in place_ids:
            plt.plot(dict_of_tokens.get(x), label=x)
        
        plt.legend(fontsize=10) #brandon
        plt.xlabel('Time-step')
        plt.ylabel('Mean tokens')
        plt.show()

def smoothen(array, filter_size):
    filt=np.ones(filter_size)/filter_size
    return convolve(array[:-(filter_size-1),0],filt)


def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

    
def create_plot(analysis, input_place_list, place_labels, mutation_list, mutation_labels, plot_title):
    

    t=np.arange(0,(desired_plotting_steps/1000),0.001) #divide middle number by 0.001 to get ur number of time steps

    t=np.arange(0,6000,0.001) #divide middle number by 0.001 to get ur number of time steps
    t=np.arange(0,(desired_plotting_steps/1000),0.001) #divide middle number by 0.001 to get ur number of time steps
    fig,ax=plt.subplots()
    linestep = 0.3
    line_width = 3
    
    for i, mutation in enumerate(mutation_list):
        for place, place_label in zip(input_place_list, place_labels):

            # data = analysis[mutation].mean_token_history_for_places([place])[0:6000000]
            # meandata = data[5000000:6000000]
            # print(sum(meandata)/len(meandata))
            # print(data[510000])
            data = analysis[mutation].mean_token_history_for_places([place])[0:desired_plotting_steps]
            y = data[:,0]
            # print(data[3000])
            if place_label == "":
                ax.plot(t, data, color = "dimgrey", label = mutation_labels[i], linewidth = line_width- i*linestep)
                y_av = movingaverage(y, 100000)
                ax.plot(t[1000:], y_av[1000:], label = 'rolling average', linewidth = line_width- i*linestep, color = "r")
                ylim(min(y), max(y))
            else:
                ax.plot(t, data, label = mutation_labels[i]+' - '+place_label, linewidth = line_width- i*linestep)
               
    
    ax.legend()
    Analysis.standardise_plot(ax, title = plot_title, xlabel = "Time (s)",ylabel = "Molecule count")
    plt.show()
    


def plot_stacked_bars(ax, legend, all_data, xs, labels, width):

    cum_sum = np.zeros_like(all_data[:,0])
    for i in range(len(labels)):
        data = all_data[:,i]
        rects = ax.bar(xs, data, width, bottom=cum_sum, label=labels[i])
        cum_sum += data    
    
def create_bar_chart(analysis, places_a, places_a_labels, places_b, places_b_labels, mutation_list, mutation_labels, plot_title):
#     for mutation in mutation_list:
#         for place in places_a:
#             print(place)
#             print(analysis[mutation].mean_token_history_for_places(place)[-1])
#     for mutation in mutation_list:
#         for place in places_b:
#             print(place)
#             print(analysis[mutation].mean_token_history_for_places(place)[-1])
    final_token_count_a = [[analysis[mutation].mean_token_history_for_places(place)[-1] for place in places_a] for mutation in mutation_list]
    final_token_count_b = [[analysis[mutation].mean_token_history_for_places(place)[-1] for place in places_b] for mutation in mutation_list]
    print(np.array(final_token_count_a).shape)
    print(np.array(final_token_count_b).shape)
    final_token_count_a = np.sum(final_token_count_a, 2) # remove dimension 3
    final_token_count_b = np.sum(final_token_count_b, 2) # remove dimension 3

    # normalize data

    final_token_count_a = final_token_count_a / np.sum(final_token_count_a[0,:])
    final_token_count_b = final_token_count_b / np.sum(final_token_count_b, 1)[:,None]

    final_token_count_a *= 100
    final_token_count_b *= 100
    
    width = 0.5
    
    FIGURESIZE = (14,7)
    fig, ax = plt.subplots(1, 1, figsize=FIGURESIZE)

    bar_positions_a = np.array(range(len(mutation_list)))
    bar_positions_b = max(bar_positions_a) + 2 + np.array(range(len(mutation_list)))
    
    plot_stacked_bars(ax,legend=mutation_list, all_data=final_token_count_a, xs=bar_positions_a, labels=places_a_labels,width=width)
    plot_stacked_bars(ax,legend=mutation_list, all_data=final_token_count_b, xs=bar_positions_b, labels = places_b_labels,width=width)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('% Molecule Count', fontsize=16)
    ax.set_title(plot_title, fontsize=18)
    ax.set_xticks(np.concatenate((bar_positions_a, bar_positions_b)))
    ax.set_xticklabels(np.concatenate((mutation_labels, mutation_labels)), rotation=-25, ha='left', fontsize=12)

    #ax.set_ylim((0,150))

    plt.legend(fontsize=14, loc='upper right', bbox_to_anchor=(1.3, 1))
    plt.show()


# # Bar charts

# In[ ]:


# create_bar_chart(analysis, 
#                  places_a = ['p_ROS_mito', 'p_Ca_cyto'], 
#                  places_a_labels = ['p_ROS_mito', 'p_Ca_cyto'], 
#                  places_b = ['p_RTN3_HMW_dys1','p_RTN3_HMW_dys2','p_RTN3_HMW_lyso'], 
#                  places_b_labels=['Dystrophic neurites I', 'Dystrophic_neurites II','Lyso'], 
#                  mutation_list = ['healthy', 'chol600'], 
#                  mutation_labels = ['healthy', 'chol600'],
#                  plot_title = 'PD - RTN3 distribution')


# In[ ]:


# create_bar_chart(analysis, 
#                  places_a = ['p_RTN3_axon','p_RTN3_PN'], 
#                  places_a_labels = ['Axon', 'Perinuclear region'], 
#                  places_b = ['p_RTN3_HMW_dys1','p_RTN3_HMW_dys2','p_RTN3_HMW_lyso'], 
#                  places_b_labels=['Dystrophic neurites I', 'Dystrophic_neurites II','Lyso'], 
#                  mutation_list = ['all_mutations', 'lrrk2', 'DNL','NPT','LAMP2A', 'healthy'], 
#                  mutation_labels = ['Combined diseased state','LRRK2','LRRK2 + DNL151','Combined diseased state + NPT200','Combined diseased state + LAMP2A', 'Healthy'],
#                  plot_title = 'PD - RTN3 distribution and therapeutics')


# # Plotting

# ## Energy metabolism 

# In[4]:


# create_plot(analysis, 
#             input_place_list = ['p_tauP'], 
#             place_labels = [""], 
#             mutation_list = ['Whole_Module_10e6'], 
#             mutation_labels = ['Whole_Module_10e6'],
#             plot_title = 'PD - p_tauP')
# create_plot(analysis, 
#             input_place_list = ['p_tauP'], 
#             place_labels = [""], 
#             mutation_list = ['Whole_Module_10e6_SD_10_percent'], 
#             mutation_labels = ['Whole_Module_10e6_SD_10_percent'],
#             plot_title = 'PD - p_tauP')

create_plot(analysis, 
            input_place_list = ['p_Ab'], 
            place_labels = [""], 
            mutation_list = ['healthy6t'], 
            mutation_labels = ['healthy'],
            plot_title = 'Ab')



#Whole_Module_10e6_SD_10_percent
#

# create_plot(analysis, 
#             input_place_list = ['p_cas3'], 
#             place_labels = [""], 
#             mutation_list = ['testlol'], 
#             mutation_labels = ['testlol'],
#             plot_title = 'PD - p_cas3')

# create_plot(analysis, 
#             input_place_list = ['p_H2O_mito'], 
#             place_labels = [""], 
#             mutation_list = ['testlol'], 
#             mutation_labels = ['testlol'],
#             plot_title = 'PD - p_H2O_mito')

# create_plot(analysis, 
#             input_place_list = ['p_chol_mito'], 
#             place_labels = [""], 
#             mutation_list = ['testlol'], 
#             mutation_labels = ['testlol'],
#             plot_title = 'PD - p_chol_mito')

# # ## Lewy body formation

# # In[ ]:


# create_plot(analysis, 
#             input_place_list = ['p_ATP'], 
#             place_labels = [""], 
#             mutation_list = ['healthy', 'chol600'], 
#             mutation_labels = ['healthy', 'chol600'],
#             plot_title = 'PD - Lewy body formation')


# # ## Chol (LB and cas3)

# # In[ ]:


# #THE CORRECT ONE FOR CHOL
# create_plot(analysis, 
#             input_place_list = ['p_LB'], 
#             place_labels = [""], 
#             mutation_list = ['healthy','gba1_lrrk2','27OHchol','27OH_lrrk2_gba1','ApoEchol','ApoE_lrrk2_gba1'], 
#             mutation_labels = ['Healthy','GBA1 + LRRK2','2x 27OH-chol','2x 27OH-chol + LRRK2 + GBA1','2x APOE-chol','2x APOE-chol + LRRK2 + GBA1'],
#             plot_title = 'PD - Lewy body formation and high levels chol')
# #THE CORRECT ONE FOR CHOL
# create_plot(analysis, 
#             input_place_list = ['p_cas3'], 
#             place_labels = [""], 
#             mutation_list = ['healthy','gba1_lrrk2','27OHchol','27OH_lrrk2_gba1','ApoEchol','ApoE_lrrk2_gba1'], 
#             mutation_labels = ['Healthy','GBA1 + LRRK2','2x 27OH-chol','2x 27OH-chol + LRRK2 + GBA1','2x APOE-chol','2x APOE-chol + LRRK2 + GBA1'],
#             plot_title = 'PD - Active Caspase-3 and high levels chol')


# # ## Therapeutics

# # In[ ]:


# create_plot(analysis, 
#             input_place_list = ['p_cas3'], 
#             place_labels = [""], 
#             mutation_list = ['all_mutations', 'lrrk2', 'DNL','NPT','LAMP2A', 'healthy'], 
#             mutation_labels = ['Combined diseased state','LRRK2','LRRK2 + DNL151','Combined diseased state + NPT200','Combined diseased state + LAMP2A', 'Healthy'],
#             plot_title = 'PD - Active Caspase-3 and therapeutics')
# # create_plot(analysis, 
# #             input_place_list = ['p_SNCA_olig'], 
# #             place_labels = [""], 
# #             mutation_list = ['all_mutations', 'lrrk2', 'DNL','NPT','LAMP2A','healthy'], 
# #             mutation_labels = ['Combined diseased state','LRRK2','LRRK2 + DNL151','Combined diseased state + NPT200','Combined diseased state + LAMP2A', 'Healthy'],
# #             plot_title = 'PD - SNCA oligomerisation and therapeutics')
# create_plot(analysis, 
#             input_place_list = ['p_LB'], 
#             place_labels = [""], 
#             mutation_list = ['all_mutations', 'lrrk2', 'DNL','NPT','LAMP2A','healthy'], 
#             mutation_labels = ['Combined diseased state','LRRK2','LRRK2 + DNL151','Combined diseased state + NPT200','Combined diseased state + LAMP2A', 'Healthy'],
#             plot_title = 'PD - Lewy body formation and therapeutics')


# # # Computing the mean

# # In[ ]:


# mean_healthy = np.mean(analysis['healthy'].token_storage[:,50000:,analysis['healthy'].place_dict["p_ATP"]])
# print("healthy", mean_healthy)
# mean_lrrk2 = np.mean(analysis['lrrk2'].token_storage[:,50000:,analysis['lrrk2'].place_dict["p_ATP"]])
# print("lrrk2", mean_lrrk2)


# # In[ ]:


# # create_plot(['p_LB'],"Lewy body formation")


# # In[ ]:


# # create_plot(['p_SNCA_olig'],"SNCA olgiomerisation")


# # In[ ]:





# # In[ ]:


# create_plot(['p_chol_LE'],"Cholesterol late endosomes")


# In[ ]:





