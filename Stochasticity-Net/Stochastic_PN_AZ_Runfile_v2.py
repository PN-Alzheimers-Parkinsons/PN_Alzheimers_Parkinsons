# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 15:42:54 2020

@author: brand
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



def main():
    
        #only runs this chunk of code if running the file directly, not as an import

    
        # Initialize an empty Petri net
    pn = PetriNet(number_of_runs=1) # so when you assign PetriNet class to pn, you also assign PetriNetModel class to the object petri_net_model since they are linked that way. This occurs because pn goes via "self" and what actually happens, is that you get pn.petri_net_model as an attribute?(right word?s)
    
#         # Add places for each chemical species

   #AB pathology
    pn.add_place(200, place_id="p_asec", label="alpha secretase")
    # pn.add_place(it_p_APP_pm, place_id="p_APP_PM", label="APP at plasma membrane")
    # pn.add_place(it_p_APP_endo, place_id="p_APP_endo", label="endocytosed APP")
    # pn.add_place(it_p_sAPPa, place_id="p_sAPPa", label="soluble sAPP alpha")
    # pn.add_place(it_p_CTF83, place_id="p_CTF83", label="CTF83")
    # pn.add_place(it_p_bsec, place_id="p_bsec", label="beta secretase")
    # pn.add_place(it_p_sAPPb, place_id="p_sAPPb", label="soluble sAPP beta")
    # pn.add_place(it_p_CTF99, place_id="p_CTF99", label="CTF99")
    # pn.add_place(it_p_Ab, place_id="p_AB", label="Amyloid beta peptide")
    # pn.add_place(it_p_AICD, place_id="p_AICD", label="AICD")
    # pn.add_place(it_p_gsec, place_id="p_gsec", label="gamma secretase")

    
    #AB pathology transitions
    # pn.add_transition(transition_id = 't_asec_exp',
    #                 label                = "alpha secretase expression",
    #                 input_place_ids       = [],
    #                 input_arc_weights  = [], 
    #                 output_place_ids       = ['p_asec'],
    #                 output_arc_weights = [1],
    #                 distribution_type = ["g", 0, 1])
   
    pn.add_transition(transition_id = 't_asec_deg',
                   label      =     "alpha secretase degradation",
                   input_place_ids         =  ['p_asec'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  [],
                   output_arc_weights =  [],
                   distribution_type = ["grf", 10, 1, r_t_asec_degr])
   
    # pn.add_transition(transition_id = 't_APP_asec_cleav',
    #                label      =     "APP cleavage by alpha secretase",
    #                input_place_ids         =  ['p_APP_PM'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_sAPPa', 'p_CTF83'],
    #                output_arc_weights =  [1, 1],
    #                catal_place_ids = ["p_asec"],
    #                catal_arc_weights = [1],
    #                distribution_type=["g",4,1,vmax_scaling_t_APP_asec_cleav]) #middle value = rate
    

    #changed catal_arc_weight to 1 so that there is catalysis occurring brandon #also changing the catal_arc_weight doesnt change very much, oh wait, might need to add a new argument called catal_arc_weight, which will change the threshold of catalysis. Right now, the threshold of catalysis is determined by the ARC WEIGHT instead of the catal ARC WEIGHT
           
    # pn.add_transition(transition_id = 't_APP_endo',
    #                label      =     "APP endocytosis",
    #                input_place_ids         =  ['p_APP_PM'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_APP_endo'],
    #                output_arc_weights =  [1])
               
    # pn.add_transition(transition_id = 't_APP_endo_bsec_cleav',
    #                label      =     "APP cleavage by beta secretase",
    #                input_place_ids         =  ['p_APP_endo'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_sAPPb', 'p_CTF99'],
    #                output_arc_weights =  [1, 1],
    #                catal_place_ids = ['p_bsec'],
    #                catal_arc_weights= [1])#here
   
   
    # pn.add_transition(transition_id = 't_bsec_exp',
    #                label      =     "beta secretase expression",
    #                input_place_ids         =  [],
    #                input_arc_weights  =  [],
    #                output_place_ids         =  ['p_bsec'],
    #                output_arc_weights =  [1])
       
    # pn.add_transition(transition_id = 't_bsec_deg',
    #                label      =     "beta secretase degradation",
    #                input_place_ids         =  ['p_bsec'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  [],
    #                output_arc_weights =  [])
       
    # pn.add_transition(transition_id = 't_CTF99_gsec_cleav',
    #                label      =     "CTF99 cleavage by gamma secretase",
    #                input_place_ids         =  ['p_CTF99'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_AB', 'p_AICD'],
    #                output_arc_weights =  [1, 1],
    #                catal_place_ids = ["p_gsec"],
    #                catal_arc_weights = [1])
       
    # pn.add_transition(transition_id = 't_gsec_exp',
    #                label      =     "gamma secretase expression",
    #                input_place_ids         =  [],
    #                input_arc_weights  =  [],
    #                output_place_ids         =  ['p_gsec'],
    #                output_arc_weights =  [1])
       
    # pn.add_transition(transition_id = 't_gsec_deg',
    #                label      =     "gamma secretase degradation",
    #                input_place_ids         =  ['p_gsec'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  [],
    #                output_arc_weights =  [])
       
    # pn.add_transition(transition_id = 't_bsec_deg',
    #                label      =     "beta secretase degradation",
    #                input_place_ids         =  ['p_bsec'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  [],
    #                output_arc_weights =  [])
    
    
    # #Tau pathology
    
    # pn.add_place(initial_tokens=0, place_id="p_GSK3b_inact", label="Inactive GSK3 beta kinase")
    # pn.add_place(1, place_id="p_GSK3b_act", label="Active GSK3 beta kinase")
    # pn.add_place(0, place_id="p_tauP", label="Phosphorylated tau")
    # pn.add_place(1, place_id="p_tau", label="Unphosphorylated tau (microtubule)")

    # #Tau pathology transitions
    
    # pn.add_transition(transition_id = 't_GSK3b_exp',
    #                label                = "GSK3 beta expression",
    #                input_place_ids       = [],
    #                input_arc_weights  = [], 
    #                output_place_ids       = ['p_GSK3b_inact'],
    #                output_arc_weights = [1])
   
    # pn.add_transition(transition_id = 't_GSK3b_deg',
    #                label      =     "GSK3 beta degradation",
    #                input_place_ids         =  ['p_GSK3b_inact'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  [],
    #                output_arc_weights =  [])
   
    # pn.add_transition(transition_id = 't_actv_GSK3b',
    #                label      =     "GSK3 beta activation",
    #                input_place_ids         =  ['p_GSK3b_inact'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_GSK3b_act'],
    #                output_arc_weights =  [1])
                 
    # pn.add_transition(transition_id = 't_inactv_GSK3b',
    #                label      =     "GSK3 beta inactivation",
    #                input_place_ids         =  ['p_GSK3b_act'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_GSK3b_inact'],
    #                output_arc_weights =  [1])
    
    # pn.add_transition(transition_id = 't_phos_tau',
    #                label      =     "Phosphorylation of tau",
    #                input_place_ids         =  ['p_tau'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_tauP'],
    #                output_arc_weights =  [1],
    #                catal_place_ids= ['p_GSK3b_act'],
    #                catal_arc_weights= [1])
    
    # pn.add_transition(transition_id = 't_dephos_tau',
    #                label      =     "Dephosphorylation of tau",
    #                input_place_ids         =  ['p_tauP'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_tau'],
    #                output_arc_weights =  [1])
    
    # pn.add_transition(transition_id = 't_inactv_GSK3b',
    #                label      =     "GSK3 beta inactivation",
    #                input_place_ids         =  ['p_GSK3b_act'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_GSK3b_inact'],
    #                output_arc_weights =  [1])
            
    
    # #Cholesterol homeostasis 
    # # Cholesterol-ApoE
    # pn.add_place(10, place_id="p_ApoEchol_extra", label="ApoE-cholesterol complex (extracellular)")
    # #Cholesterol in different organelles
    # pn.add_place(10, place_id="p_chol_LE", label="Cholesterol (late endosome)")
    # pn.add_place(1, place_id="p_chol_mito", label="Cholesterol (mitochondria)")
    # pn.add_place(1, place_id="p_chol_ER", label="Cholesterol (ER)")
    # pn.add_place(1, place_id="p_chol_PM", label="Cholesterol (plasma membrane)") 
    #  # Oxysterols
    # pn.add_place(1, place_id="p_24OHchol_extra", label="24-hydroxycholesterol (extracellular)")
      # pn.add_place(1, place_id="p_24OHchol_intra", label="24-hydroxycholesterol (intracellular)")
    # pn.add_place(1, place_id="p_27OHchol_extra", label="27-hydroxycholesterol (extracellular)")
    # pn.add_place(1, place_id="p_27OHchol_intra", label="27-hydroxycholesterol (intracellular)")
    # pn.add_place(1, place_id="p_7HOCA", label="7-HOCA")
    # pn.add_place(1, place_id="p_preg", label="Pregnenolon")
    
    # #Transitions
    # pn.add_transition(transition_id = 't_LDLR_endocyto',
    #                label      =     "LDLR endocytosis",
    #                input_place_ids         =  ['p_ApoEchol_extra'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_chol_LE'],
    #                output_arc_weights =  [1],
    #                inhib_place_ids= ['p_chol_ER'],
    #                inhib_arc_weights= [1]) #not sure what to use for this 
    
    # pn.add_transition(transition_id = 't_chol_trans_LE_ER',
    #                label      =     "Chol transport LE-ER",
    #                input_place_ids         =  ['p_chol_LE'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_chol_ER'],
    #                output_arc_weights =  [1])
    
    # pn.add_transition(transition_id = 't_chol_trans_LE_mito',
    #                label      =     "Chol transport LE-mito",
    #                input_place_ids         =  ['p_chol_LE'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_chol_mito'],
    #                output_arc_weights =  [1])
    
    # pn.add_transition(transition_id = 't_chol_trans_LE_PM',
    #                label      =     "Chol transport LE-PM",
    #                input_place_ids         =  ['p_chol_LE'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_chol_PM'],
    #                output_arc_weights =  [1])
    
    # pn.add_transition(transition_id = 't_chol_trans_PM_ER',
    #                label      =     "Chol transport PM-ER",
    #                input_place_ids         =  ['p_chol_PM'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_chol_ER'],
    #                output_arc_weights =  [1])
    
    # pn.add_transition(transition_id = 't_chol_trans_ER_mito',
    #                label      =     "Chol transport ER-mito",
    #                input_place_ids         =  ['p_chol_ER'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_chol_mito'],
    #                output_arc_weights =  [1])
    
    # pn.add_transition(transition_id = 't_CYP27A1_metab',
    #                label      =     "Chol metab CYP27A1",
    #                input_place_ids         =  ['p_chol_mito'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_27OHchol_intra'],
    #                output_arc_weights =  [1])                                          
 
    # pn.add_transition(transition_id = 't_CYP11A1_metab',
    #                label      =     "Chol metab CYP11A1",
    #                input_place_ids         =  ['p_chol_mito'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_preg'],
    #                output_arc_weights =  [1])  
                                                
    # pn.add_transition(transition_id = 't_CYP7B1_metab',
    #                label      =     "27OHchol metab CYP7B1",
    #                input_place_ids         =  ['p_27OHchol_intra'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_7HOCA'],
    #                output_arc_weights =  [1])                                                  
    
    # pn.add_transition(transition_id = 't_27OHchol_endocytp',
    #                label      =     "27OHchol endocyto",
    #                input_place_ids         =  ['p_27OHchol_extra'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_27OHchol_extra', 'p_27OHchol_intra'],
    #                output_arc_weights =  [1,1])  
                                                
    # pn.add_transition(transition_id = 't_CYP46A1_metab',
    #                label      =     "Chol metab CYP46A1",
    #                input_place_ids         =  ['p_chol_ER'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_24OHchol_intra'],
    #                output_arc_weights =  [1])                                                  
     
    # pn.add_transition(transition_id = 't_24OHchol_exocyto',
    #                label      =     "24OHchol exocyto",
    #                input_place_ids         =  ['p_24OHchol_intra'],
    #                input_arc_weights  =  [1],
    #                output_place_ids         =  ['p_24OHchol_extra'],
    #                output_arc_weights =  [1])                                                  
      
    # pn.add_transition(transition_id = 't_chol_trans_PM_ECM',
    #                label      =     "Chol transport PM-ECM",
    #                input_place_ids         =  ['p_chol_PM', 'p_24OHchol_intra'],
    #                input_arc_weights  =  [1, 0],
    #                output_place_ids         =  [],
    #                output_arc_weights =  [])                                                  
      
       
       # Run the network X times
    #a = {place.place_id:place.tokens for place in petri_net_model.places.values()}
    pn.run(2000, print_stats=False)
    
    #BSL: A good looking curve is, 2000 run steps, standard deviation of fixed 1 token, 200 starting tokens for asec. 10% of the mean gives a very smooth curve.

    # Plot the time-evolution of the system
    #input the place ids into this list for plotting
    list_for_plot = ['p_asec'] 
    
    pn.plot_time_evolution(list_for_plot)


    # Generate block diagram of the Petri net

    #pn.generate_diagram("/Users/brand/Documents/Github/PN_Alzheimers_Parkinsons/blockdiags/AD")

if __name__ == "__main__":
    main()
