#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 10:46:14 2021

@author: madeleinericher
"""
import os
import sys

cwd = os.getcwd() # Get current working directory
root_folder = os.sep + "PN_Alzheimers_Parkinsons"
# Move to 'utils' from current directory position
sys.path.insert(0, cwd[:(cwd.index(root_folder)+len(root_folder))] + os.sep + "Stochasticity-Net" + os.sep + "ADutils" + os.sep)

from Stochastic_PN_Architecture_v2 import *

from parameters import *
from rate_functions import *
from initial_tokens import *
from firing_conditions import *
from Stochastic_Analysis import Analysis



def main():
    
        #only runs this chunk of code if running the file directly, not as an import

    
        # Initialize an empty Petri net
    pn = PetriNet(number_of_runs=1) # so when you assign PetriNet class to pn, you also assign PetriNetModel class to the object petri_net_model since they are linked that way. This occurs because pn goes via "self" and what actually happens, is that you get pn.petri_net_model as an attribute?(right word?s)
    
#         # Add places for each chemical species

   #AB pathology
    pn.add_place(it_p_asec, place_id="p_asec", label="alpha secretase")
    pn.add_place(it_p_APP_pm, place_id="p_APP_pm", label="APP at plasma membrane")
    pn.add_place(it_p_APP_endo, place_id="p_APP_endo", label="endocytosed APP")
    pn.add_place(it_p_sAPPa, place_id="p_sAPPa", label="soluble sAPP alpha")
    pn.add_place(it_p_CTF83, place_id="p_CTF83", label="CTF83")
    pn.add_place(it_p_bsec, place_id="p_bsec", label="beta secretase")
    pn.add_place(it_p_sAPPb, place_id="p_sAPPb", label="soluble sAPP beta")
    pn.add_place(it_p_CTF99, place_id="p_CTF99", label="CTF99")
    pn.add_place(it_p_Ab, place_id="p_Ab", label="Amyloid beta peptide")
    pn.add_place(it_p_AICD, place_id="p_AICD", label="AICD")
    pn.add_place(it_p_gsec, place_id="p_gsec", label="gamma secretase")
    pn.add_place(it_p_age, place_id="p_age", label="age risk factor")
    pn.add_place(it_p_ApoE, place_id="p_ApoE", label="ApoE risk factor")
  
    pn.add_place(it_p_chol_LE, place_id="p_chol_LE", label="Cholesterol (late endosome)")
    pn.add_place(it_p_chol_mito, place_id="p_chol_mito", label="Cholesterol (mitochondria)")
    pn.add_place(it_p_chol_ER, place_id="p_chol_ER", label="Cholesterol (ER)")
    pn.add_place(it_p_chol_PM, place_id="p_chol_PM", label="Cholesterol (plasma membrane)") 
    #   # Oxysterols
    pn.add_place(it_p_24OHchol_extra, place_id="p_24OHchol_extra", label="24-hydroxycholesterol (extracellular)")
    pn.add_place(it_p_24OHchol_intra, place_id="p_24OHchol_intra", label="24-hydroxycholesterol (intracellular)")
    pn.add_place(it_p_27OHchol_extra, place_id="p_27OHchol_extra", label="27-hydroxycholesterol (extracellular)")
    pn.add_place(it_p_27OHchol_intra, place_id="p_27OHchol_intra", label="27-hydroxycholesterol (intracellular)")
    pn.add_place(it_p_cas3, place_id="p_cas3", label="Active caspase 3")
    pn.add_place(it_p_ATP, place_id="p_ATP", label="ATP")
    pn.add_place(it_p_ADP, place_id="p_ADP", label="ADP")
    pn.add_place(it_p_reduc_mito, place_id="p_reduc_mito", label="Reducing agents (mitochondria)")
    pn.add_place(it_p_ROS_mito, place_id="p_ROS_mito", label="ROS (mitochondria)")
    pn.add_place(it_p_H2O_mito, place_id="p_H2O_mito", label="H2O (mitochondria)")
    pn.add_place(it_p_RTN3_axon, place_id="p_RTN3_axon", label="Monomeric RTN3 (axonal)")
    pn.add_place(it_p_RTN3_PN, place_id="p_RTN3_PN", label="Monomeric RTN3 (perinuclear)")

    # HMW RTN3 (cycling between different cellular compartments)
    pn.add_place(it_p_RTN3_HMW_cyto, place_id="p_RTN3_HMW_cyto", label="HMW RTN3 (cytosol)")
    pn.add_place(it_p_RTN3_HMW_auto, place_id="p_RTN3_HMW_auto", label="HMW RTN3 (autophagosome)")
    pn.add_place(it_p_RTN3_HMW_lyso, place_id="p_RTN3_HMW_lyso", label="HMW RTN3 (degraded in lysosome)")
    pn.add_place(it_p_RTN3_HMW_dys1, place_id="p_RTN3_HMW_dys1", label="HMW RTN3 (type I/III dystrophic neurites)")
    pn.add_place(it_p_RTN3_HMW_dys2, place_id="p_RTN3_HMW_dys2", label="HMW RTN3 (type II dystrophic neurites)")
    
    
    ##AB aggregation places
    pn.add_place(it_p_Ab_elon, place_id="p_Ab_elon", label="Elongating Ab")
    pn.add_place(it_p_Ab_olig, place_id="p_Ab_olig", label="Ab oligomer")
    pn.add_place(it_p_Ab_fib, place_id="p_Ab_fib", label="Ab fibril")

    
    # AB pathology transitions
    pn.add_transition(transition_id = 't_asec_exp',
                        label                = "alpha secretase expression",
                        input_place_ids       = ['p_24OHchol_intra'],
                        input_arc_weights  = [0], 
                        output_place_ids       = ['p_asec'],
                        output_arc_weights = [1],
                        distribution_type = ["grf", SD, r_t_asec_exp, fc_t_asec_exp])
       
    pn.add_transition(transition_id = 't_asec_deg',
                        label      =     "alpha secretase degradation",
                        input_place_ids         =  ['p_asec'],
                        input_arc_weights  =  [1],
                        output_place_ids         =  [],
                        output_arc_weights =  [],
                        distribution_type = ["grf", SD, r_t_asec_degr, fc_t_asec_degr])
   
    pn.add_transition(transition_id = 't_APP_asec_cleav',
                        label      =     "APP cleavage by alpha secretase",
                        input_place_ids         =  ['p_APP_pm', 'p_asec', 'p_chol_PM'],
                        input_arc_weights  =  [1, 0, 0],
                        output_place_ids         =  ['p_sAPPa', 'p_CTF83'],
                        output_arc_weights =  [1, 1],
                        distribution_type=["grf", SD, r_t_APP_asec_cleav, fc_t_APP_asec_cleav]) 
    
    pn.add_transition(transition_id                 = 't_APP_exp',
                        label                         = 'APP expression rate',
                        input_place_ids                 = ['p_ApoE', 'p_ROS_mito'],
                        input_arc_weights     = [0, 0], 
                        output_place_ids = ['p_APP_pm'],
                        output_arc_weights = [1],
                        distribution_type = ["grf", SD, r_t_APP_exp, fc_t_APP_exp]) 
            
    pn.add_transition(transition_id                 = 't_APP_endo_event',
                        label                         = 'APP-utilizing cellular events',
                        input_place_ids                 = ['p_APP_endo'], 
                        input_arc_weights     = [1], 
                        output_place_ids = [],
                        output_arc_weights = [],
                        distribution_type = ["grf", SD, r_t_APP_endo_event, fc_t_APP_endo_event])
    
    pn.add_transition(transition_id = 't_APP_endo',
                        label      =     "APP endocytosis",
                        input_place_ids         =  ['p_APP_pm', 'p_ApoE'],
                        input_arc_weights  =  [1, 0],
                        output_place_ids         =  ['p_APP_endo'],
                        output_arc_weights =  [1],
                        distribution_type = ["grf", SD, r_t_APP_endocyto, fc_t_APP_endocyto])
                   
    pn.add_transition(transition_id = 't_APP_bsec_cleav',
                        label      =     "APP cleavage by beta secretase",
                        input_place_ids         =  ['p_APP_endo', 'p_bsec', 'p_chol_PM', 'p_age'],
                        input_arc_weights  =  [1, 0, 0, 0],
                        output_place_ids         =  ['p_sAPPb', 'p_CTF99'],
                        output_arc_weights =  [1, 1],
                        distribution_type = ["grf", SD, r_t_APP_bsec_cleav, fc_t_APP_bsec_cleav])
   
    pn.add_transition(transition_id = 't_bsec_exp',
                        label      =     "beta secretase expression",
                        input_place_ids         =  ['p_ROS_mito', 'p_27OHchol_intra', 'p_RTN3_axon'],
                        input_arc_weights  =  [0, 0, 0],
                        output_place_ids         =  ['p_bsec'],
                        output_arc_weights =  [1],
                        distribution_type = ["grf", SD, r_t_bsec_exp, fc_t_bsec_exp])
       
    pn.add_transition(transition_id = 't_bsec_deg',
                        label      =     "beta secretase degradation",
                        input_place_ids         =  ['p_bsec'],
                        input_arc_weights  =  [1],
                        output_place_ids         =  [],
                        output_arc_weights =  [],
                        distribution_type = ["grf", SD, r_t_bsec_degr, fc_t_bsec_degr])
       
    pn.add_transition(transition_id = 't_CTF99_gsec_cleav',
                        label      =     "CTF99 cleavage by gamma secretase",
                        input_place_ids         =  ['p_CTF99', 'p_chol_PM'],
                        input_arc_weights  =  [1, 0],
                        output_place_ids         =  ['p_Ab', 'p_AICD'],
                        output_arc_weights =  [1, 1],
                        catal_place_ids = ["p_gsec"],
                        catal_arc_weights = [1],
                        distribution_type = ["grf", SD, r_t_CTF99_gsec_cleav, fc_CTF99_gsec_cleav])
       
    pn.add_transition(transition_id = 't_gsec_exp',
                        label      =     "gamma secretase expression",
                        input_place_ids         =  ['p_ROS_mito'],
                        input_arc_weights  =  [0],
                        output_place_ids         =  ['p_gsec'],
                        output_arc_weights =  [1],
                        distribution_type = ["grf", SD, r_t_gsec_exp, fc_t_gsec_exp])
           
    pn.add_transition(transition_id = 't_gsec_deg',
                        label      =     "gamma secretase degradation",
                        input_place_ids         =  ['p_gsec'],
                        input_arc_weights  =  [1],
                        output_place_ids         =  [],
                        output_arc_weights =  [],
                        distribution_type = ["grf", SD, r_t_gsec_degr, fc_t_gsec_degr])

    pn.add_transition(transition_id                 = 't_Ab_degr',
                        label                         = 'Ab degradation',
                        input_place_ids                 = ['p_Ab'], 
                        input_arc_weights     = [1], 
                        output_place_ids = [],
                        output_arc_weights = [],
                        distribution_type = ["grf", SD, r_t_Ab_degr, fc_t_Ab_degr]) # TODO - fix ratio    
    
    
    
    ##AB Aggregation transitions
    
    pn.add_transition(transition_id = 't_Ab_elon',
                        label                = "Ab elongation step",
                        input_place_ids       = ['p_Ab'],
                        input_arc_weights  = [1], 
                        output_place_ids       = ['p_Ab_elon'],
                        output_arc_weights = [1],
                        distribution_type = ["grf", 10, r_t_Ab_elon, fc_t_Ab_elon])

    pn.add_transition(transition_id = 't_Ab_agg',
                        label                = "Ab aggregation",
                        input_place_ids       = ['p_Ab_elon'],
                        input_arc_weights  = [12.4], 
                        output_place_ids       = ['p_Ab_olig'],
                        output_arc_weights = [1],
                        distribution_type = ["grf", 0, r_t_Ab_agg, fc_t_Ab_agg])

    
    pn.add_transition(transition_id = 't_Ab_fib',
                        label                = "Ab fibrillation",
                        input_place_ids       = ['p_Ab_olig'],
                        input_arc_weights  = [4], 
                        output_place_ids       = ['p_Ab_fib'],
                        output_arc_weights = [1],
                        distribution_type = ["grf", 0, r_t_Ab_fib, fc_t_Ab_fib])

    
    
    # #Simpler model
    # pn.add_transition(transition_id = 't_Ab_agg3',
    #                     label                = "Ab aggregation",
    #                     input_place_ids       = ['p_Ab'],
    #                     input_arc_weights  = [12.4], 
    #                     output_place_ids       = ['p_Ab_olig'],
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0, r_t_Ab_agg3, fc_t_Ab_elon])

    # pn.add_transition(transition_id = 't_Ab_fib3',
    #                     label                = "Ab fibrillation",
    #                     input_place_ids       = ['p_Ab_olig'],
    #                     input_arc_weights  = [4], 
    #                     output_place_ids       = ['p_Ab_fib'],
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0, r_t_Ab_agg3, fc_t_Ab_fib])
    
    
    
       # Run the network X times
    #a = {place.place_id:place.tokens for place in petri_net_model.places.values()}
    pn.run(100000, print_stats=False)
    

    # Plot the time-evolution of the system
    #input the place ids into this list for plotting
    list_for_plot = ['p_Ab_elon'] 
    
    pn.plot_time_evolution(list_for_plot)
    # pn.timeseries_mean_for_place("p_Ca_extra")
    analysis = Analysis(pn)
    run_save_name = "test"
    Analysis.store_to_file(analysis, run_save_name)
    print('Network saved to : "' + run_save_name+'.pkl"')

    
if __name__ == "__main__":
    main()