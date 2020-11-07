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
sys.path.insert(0, cwd[:(cwd.index(root_folder)+len(root_folder))] + os.sep + "Stochasticity-Net" + os.sep)

from Stochastic_PN_Architecture import *

def main():
    
        #only runs this chunk of code if running the file directly, not as an import

    
        # Initialize an empty Petri net
    pn = PetriNet(number_of_runs=5) # so when you assign PetriNet class to pn, you also assign PetriNetModel class to the object petri_net_model since they are linked that way. This occurs because pn goes via "self" and what actually happens, is that you get pn.petri_net_model as an attribute?(right word?s)
    
#         # Add places for each chemical species




    # Cholesterol homeostasis
    pn.add_place(0, place_id="p_chol_PM",label="Chol - perinuclear region")
    pn.add_place(0, place_id="p_chol_LE",label="Chol - late endosome")
    pn.add_place(0, place_id="p_chol_ER",label="Chol - ER")
    pn.add_place(0, place_id="p_chol_mito",label="Chol - mitochondria")
    pn.add_place(0, place_id="p_27OHchol_extra",label="27-OH chol - extracellular")
    pn.add_place(0, place_id="p_27OHchol_intra",label="27-OH chol - intracellular")
    pn.add_place(0, place_id="p_ApoEchol_extra",label="ApoE - extracellular")
    pn.add_place(100, place_id="p_ApoEchol_EE",label="ApoE - Early endosome")
    pn.add_place(0, place_id="p_7HOCA",label="7-HOCA")
    pn.add_place(0, place_id="p_preg",label="Pregnenolon")
    pn.add_place(0, place_id="p_24OHchol_extra",label="24OHchol extra")
    pn.add_place(0, place_id="p_24OHchol_intra",label="24OHchol intra")
    
    #Lewy Body Pathology
    pn.add_place(0, place_id="p_LB", label="Lewy Body")
  
    
    #Cholesterol Homeostasis Transitions
    
    # Cholesterol Endocytosis
    pn.add_transition(transition_id = 't_LDLR_endocyto',
                   label      =     "LDLR endocyto",
                   input_place_ids         =  ['p_ApoEchol_extra', 'p_chol_ER', 'p_LB'],
                   input_arc_weights  =  [0,0,0],
                   output_place_ids         =  ['p_ApoEchol_EE'],
                   output_arc_weights =  [1],
                   distribution_type = ["g",10,3])
    
    # Cleavage of cholesteryl esters 
    pn.add_transition(transition_id = 't_ApoEchol_cleav',
                   label      =     "ApoE-chol cleav",
                   input_place_ids         =  ['p_ApoEchol_EE'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_LE'],
                   output_arc_weights =  [1], #was 354
                   distribution_type = ["g",10,3])
    # Transport Cholesterol from LE to ER
    pn.add_transition(transition_id = 't_chol_trans_LE_ER',
                   label      =     "Chol transport LE-ER",
                   input_place_ids         =  ['p_chol_LE'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_ER'],
                   output_arc_weights =  [1],
                   distribution_type = ["g",10,3])
    
    # Transport Cholesterol from LE to mito
    pn.add_transition(transition_id = 't_chol_trans_LE_mito',
                   label      =     "Chol transport LE-mito",
                   input_place_ids         =  ['p_chol_LE'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_mito'],
                   output_arc_weights =  [1],
                   distribution_type = ["g",10,3])

    #Transport Cholesterol from LE to PM
    pn.add_transition(transition_id = 't_chol_trans_LE_PM',
                    label      =     "Chol transport LE-PM",
                    input_place_ids         =  ['p_chol_LE'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_chol_PM'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])
    
    
    # Transport Cholesterol from PM to ER
    pn.add_transition(transition_id = 't_chol_trans_PM_ER',
                    label      =     "Chol transport PM-ER",
                    input_place_ids         =  ['p_chol_PM'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_chol_ER'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])    
    


    # Transport Cholesterol from ER to PM
    pn.add_transition(transition_id = 't_chol_trans_ER_PM',
                    label      =     "Chol transport ER-PM",
                    input_place_ids         =  ['p_chol_ER'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_chol_PM'],
                    output_arc_weights =  [1],
                    distribution_type = ["g",10,3])   
    


    # Transport Cholesterol from ER to mito
    pn.add_transition(transition_id = 't_chol_trans_ER_mito',
                    label      =     "Chol transport ER-mito",
                    input_place_ids         =  ['p_chol_ER'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_chol_mito'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])
    

    # Metabolisation of chol by CYP27A1
    pn.add_transition(transition_id = 't_CYP27A1_metab',
                    label      =     "Chol metab CYP27A1",
                    input_place_ids         =  ['p_chol_mito'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_27OHchol_intra'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])                                          
    
    
    # Metabolism of chol by CYP11A1
    pn.add_transition(transition_id = 't_CYP11A1_metab',
                    label      =     "Chol metab CYP11A1",
                    input_place_ids         =  ['p_chol_mito'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_preg'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])  

    # Metabolisation of 27OHchol by CYP7B1
    pn.add_transition(transition_id = 't_CYP7B1_metab',
                    label      =     "27OHchol metab CYP7B1",
                    input_place_ids         =  ['p_27OHchol_intra'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_7HOCA'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])                                                  
    
    


    # Endocytosis of 27OHchol
    pn.add_transition(transition_id = 't_27OHchol_endocyto',
                    label      =     "27OHchol endocyto",
                    input_place_ids         =  ['p_27OHchol_extra'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_27OHchol_extra', 'p_27OHchol_intra'],
                    output_arc_weights =  [1,1], 
                    distribution_type = ["g",10,3])  
    


    # Metabolisation of chol by CYP46A1
    pn.add_transition(transition_id = 't_CYP46A1_metab',
                    label      =     "Chol metab CYP46A1",
                    input_place_ids         =  ['p_chol_ER'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_24OHchol_intra'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])                                                  
    

    # Exocytosis of 24OHchol
    pn.add_transition(transition_id = 't_24OHchol_exocyto',
                    label      =     "24OHchol exocyto",
                    input_place_ids         =  ['p_24OHchol_intra'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_24OHchol_extra'],
                    output_arc_weights =  [1], 
                    distribution_type = ["g",10,3])                                                  
    

    # Transport of Chol into ECM
    pn.add_transition(transition_id = 't_chol_trans_PM_ECM',
                    label      =     "Chol transport PM-ECM",
                    input_place_ids         =  ['p_chol_PM', 'p_24OHchol_intra'],
                    input_arc_weights  =  [1, 0],
                    output_place_ids         =  [],
                    output_arc_weights =  [], 
                    distribution_type = ["g",10,3])                                                  
    





       
       # Run the network X times
    pn.run(100, print_stats=False)

    # Plot the time-evolution of the system
    #input the place ids into this list for plotting
    list_for_plot = ["p_chol_PM", "p_chol_LE", "p_chol_ER", "p_chol_mito", "p_27OHchol_extra", "p_27OHchol_intra", "p_ApoEchol_extra", "p_ApoEchol_EE","p_7HOCA", "p_preg", "p_24OHchol_extra", "p_24OHchol_intra"] 
    
    pn.plot_time_evolution(list_for_plot)


    # Generate block diagram of the Petri net

    pn.generate_diagram("/Users/brand/Documents/Github/PN_Alzheimers_Parkinsons/blockdiags/Brandon/PD")

if __name__ == "__main__":
    main()
