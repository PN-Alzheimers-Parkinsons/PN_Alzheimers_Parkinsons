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
    
    #
    pn.add_place(it_p_age, place_id="p_age", label="age risk factor")
    pn.add_place(it_p_ApoE, place_id="p_ApoE", label="ApoE risk factor")
    
    #Tau Pathology places
    
    pn.add_place(it_p_GSK3b_inact, place_id="p_GSK3b_inact", label="Inactive GSK3 beta kinase")
    pn.add_place(it_p_GSK3b_act, place_id="p_GSK3b_act", label="Active GSK3 beta kinase")
    pn.add_place(it_p_tauP, place_id="p_tauP", label="Phosphorylated tau")
    pn.add_place(it_p_tau, place_id="p_tau", label="Unphosphorylated tau (microtubule)")
  
    
    #Cholesterol homeostasis 
    # Cholesterol-ApoE
    pn.add_place(it_p_ApoEchol_extra, place_id="p_ApoEchol_extra", label="ApoE-cholesterol complex (extracellular)")
    # #Cholesterol in different organelles
    pn.add_place(it_p_chol_LE, place_id="p_chol_LE", label="Cholesterol (late endosome)")
    pn.add_place(it_p_chol_mito, place_id="p_chol_mito", label="Cholesterol (mitochondria)")
    pn.add_place(it_p_chol_ER, place_id="p_chol_ER", label="Cholesterol (ER)")
    pn.add_place(it_p_chol_PM, place_id="p_chol_PM", label="Cholesterol (plasma membrane)") 
    #   # Oxysterols
    pn.add_place(it_p_24OHchol_extra, place_id="p_24OHchol_extra", label="24-hydroxycholesterol (extracellular)")
    pn.add_place(it_p_24OHchol_intra, place_id="p_24OHchol_intra", label="24-hydroxycholesterol (intracellular)")
    pn.add_place(it_p_27OHchol_extra, place_id="p_27OHchol_extra", label="27-hydroxycholesterol (extracellular)")
    pn.add_place(it_p_27OHchol_intra, place_id="p_27OHchol_intra", label="27-hydroxycholesterol (intracellular)")
    pn.add_place(it_p_7HOCA, place_id="p_7HOCA", label="7-HOCA")
    pn.add_place(it_p_preg, place_id="p_preg", label="Pregnenolon")
    
    #ER retraction and collapse places
    pn.add_place(it_p_RTN3_axon, place_id="p_RTN3_axon", label="Monomeric RTN3 (axonal)")
    pn.add_place(it_p_RTN3_PN, place_id="p_RTN3_PN", label="Monomeric RTN3 (perinuclear)")

    # HMW RTN3 (cycling between different cellular compartments)
    pn.add_place(it_p_RTN3_HMW_cyto, place_id="p_RTN3_HMW_cyto", label="HMW RTN3 (cytosol)")
    pn.add_place(it_p_RTN3_HMW_auto, place_id="p_RTN3_HMW_auto", label="HMW RTN3 (autophagosome)")
    pn.add_place(it_p_RTN3_HMW_lyso, place_id="p_RTN3_HMW_lyso", label="HMW RTN3 (degraded in lysosome)")
    pn.add_place(it_p_RTN3_HMW_dys1, place_id="p_RTN3_HMW_dys1", label="HMW RTN3 (type I/III dystrophic neurites)")
    pn.add_place(it_p_RTN3_HMW_dys2, place_id="p_RTN3_HMW_dys2", label="HMW RTN3 (type II dystrophic neurites)")

    
    # Energy metabolism places
    pn.add_place(it_p_cas3, place_id="p_cas3", label="Active caspase 3")
    pn.add_place(it_p_ATP, place_id="p_ATP", label="ATP")
    pn.add_place(it_p_ADP, place_id="p_ADP", label="ADP")
    pn.add_place(it_p_reduc_mito, place_id="p_reduc_mito", label="Reducing agents (mitochondria)")
    pn.add_place(it_p_ROS_mito, place_id="p_ROS_mito", label="ROS (mitochondria)")
    pn.add_place(it_p_H2O_mito, place_id="p_H2O_mito", label="H2O (mitochondria)")

    # Calcium homeostasis places
    pn.add_place(it_p_Ca_cyto, place_id="p_Ca_cyto", label="Calcium (cytosol)")
    pn.add_place(it_p_Ca_mito, place_id="p_Ca_mito", label="Calcium (mitochondria)")
    pn.add_place(it_p_Ca_ER, place_id="p_Ca_ER", label="Calcium (ER)")

    # Discrete on/of-switches calcium pacemaking
    pn.add_place(1, place_id="p_Ca_extra", label="on1 - Calcium (extracellular)")
    pn.add_place(0, place_id="p_on2", label="on2")
    pn.add_place(0, place_id="p_on3", label="on3")
    pn.add_place(0, place_id="p_on4", label="on4")

    
    # AB pathology transitions
    pn.add_transition(transition_id = 't_asec_exp',
                    label                = "alpha secretase expression",
                    input_place_ids       = ['p_24OHchol_intra'],
                    input_arc_weights  = [0], 
                    output_place_ids       = ['p_asec'],
                    output_arc_weights = [1],
                    distribution_type = ["grf", SD, r_t_asec_exp])
   
    pn.add_transition(transition_id = 't_asec_deg',
                    label      =     "alpha secretase degradation",
                    input_place_ids         =  ['p_asec'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  [],
                    output_arc_weights =  [],
                    distribution_type = ["grf", SD, r_t_asec_degr]) #0.1 means 10%, 1 means 100% Standard deviation
    #relatieve standard deviation (10%, distribution and standard) #add michaeles menten
   
    pn.add_transition(transition_id = 't_APP_asec_cleav',
                    label      =     "APP cleavage by alpha secretase",
                    input_place_ids         =  ['p_APP_pm', 'p_asec', 'p_chol_PM'],
                    input_arc_weights  =  [1, 0, 0],
                    output_place_ids         =  ['p_sAPPa', 'p_CTF83'],
                    output_arc_weights =  [1, 1],
                    distribution_type=["grf", SD, r_t_APP_asec_cleav]) #middle value = rate
    

    # #changed catal_arc_weight to 1 so that there is catalysis occurring brandon #also changing the catal_arc_weight doesnt change very much, oh wait, might need to add a new argument called catal_arc_weight, which will change the threshold of catalysis. Right now, the threshold of catalysis is determined by the ARC WEIGHT instead of the catal ARC WEIGHT
    
    pn.add_transition(transition_id				 = 't_APP_exp',
                    label						 = 'APP expression rate',
                    input_place_ids				 = ['p_ApoE', 'p_ROS_mito'],
                    input_arc_weights	 = [0, 0], 
                    output_place_ids = ['p_APP_pm'],
                    output_arc_weights = [1],
                    distribution_type = ["grf", SD, r_t_APP_exp]) 
        
    pn.add_transition(transition_id				 = 't_APP_endo_event',
                    label						 = 'APP-utilizing cellular events',
                    input_place_ids				 = ['p_APP_endo'], 
                    input_arc_weights	 = [1], 
                    output_place_ids = [],
                    output_arc_weights = [],
                    distribution_type = ["grf", SD, r_t_APP_endo_event])
    
           
    pn.add_transition(transition_id = 't_APP_endo',
                    label      =     "APP endocytosis",
                    input_place_ids         =  ['p_APP_pm', 'p_ApoE'],
                    input_arc_weights  =  [1, 0],
                    output_place_ids         =  ['p_APP_endo'],
                    output_arc_weights =  [1],
                    distribution_type = ["grf", SD, r_t_APP_endocyto])
               
    pn.add_transition(transition_id = 't_APP_bsec_cleav',
                    label      =     "APP cleavage by beta secretase",
                    input_place_ids         =  ['p_APP_endo', 'p_bsec', 'p_chol_PM', 'p_age'],
                    input_arc_weights  =  [1, 0, 0, 0],
                    output_place_ids         =  ['p_sAPPb', 'p_CTF99'],
                    output_arc_weights =  [1, 1],
                    distribution_type = ["grf", SD, r_t_APP_bsec_cleav])
   
   
    pn.add_transition(transition_id = 't_bsec_exp',
                    label      =     "beta secretase expression",
                    input_place_ids         =  ['p_ROS_mito', 'p_27OHchol_intra', 'p_RTN3_axon'],
                    input_arc_weights  =  [0, 0, 0],
                    output_place_ids         =  ['p_bsec'],
                    output_arc_weights =  [1],
                    distribution_type = ["grf", SD, r_t_bsec_exp])
       
    pn.add_transition(transition_id = 't_bsec_deg',
                    label      =     "beta secretase degradation",
                    input_place_ids         =  ['p_bsec'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  [],
                    output_arc_weights =  [],
                    distribution_type = ["grf", SD, r_t_bsec_degr])
       
    pn.add_transition(transition_id = 't_CTF99_gsec_cleav',
                    label      =     "CTF99 cleavage by gamma secretase",
                    input_place_ids         =  ['p_CTF99', 'p_chol_PM'],
                    input_arc_weights  =  [1, 0],
                    output_place_ids         =  ['p_Ab', 'p_AICD'],
                    output_arc_weights =  [1, 1],
                    catal_place_ids = ["p_gsec"],
                    catal_arc_weights = [1],
                    distribution_type = ["grf", SD, r_t_CTF99_gsec_cleav])
       
    pn.add_transition(transition_id = 't_gsec_exp',
                    label      =     "gamma secretase expression",
                    input_place_ids         =  ['p_ROS_mito'],
                    input_arc_weights  =  [0],
                    output_place_ids         =  ['p_gsec'],
                    output_arc_weights =  [1],
                    distribution_type = ["grf", SD, r_t_gsec_exp])
       
    pn.add_transition(transition_id = 't_gsec_deg',
                    label      =     "gamma secretase degradation",
                    input_place_ids         =  ['p_gsec'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  [],
                    output_arc_weights =  [],
                    distribution_type = ["grf", SD, r_t_gsec_degr])

    pn.add_transition(transition_id				 = 't_Ab_degr',
                    label						 = 'Ab degradation',
                    input_place_ids				 = ['p_Ab'], 
                    input_arc_weights	 = [1], 
                    output_place_ids = [],
                    output_arc_weights = [],
                    distribution_type = ["grf", SD, r_t_Ab_degr]) # TODO - fix ratio    

    #Tau pathology

    #Tau pathology transitions
    
    pn.add_transition(transition_id = 't_GSK3b_exp_degr',
                    label                = "GSK3 beta expression and degradation",
                    input_place_ids       = ['p_GSK3b_inact'],
                    input_arc_weights  = [0], 
                    output_place_ids       = ['p_GSK3b_inact'],
                    output_arc_weights = [1],
                    distribution_type = ["grf", SD, r_t_GSK3b_exp_deg])
   
    pn.add_transition(transition_id = 't_actv_GSK3b',
                    label      =     "GSK3 beta activation",
                    input_place_ids         =  ['p_GSK3b_inact', 'p_ApoE', 'p_Ab'],
                    input_arc_weights  =  [1, 0, 0],
                    output_place_ids         =  ['p_GSK3b_act'],
                    output_arc_weights =  [1],
                    distribution_type = ["grf", SD, r_t_actv_GSK3b])
                 
    pn.add_transition(transition_id = 't_inactv_GSK3b',
                    label      =     "GSK3 beta inactivation",
                    input_place_ids         =  ['p_GSK3b_act'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_GSK3b_inact'],
                    output_arc_weights =  [1],
                    distribution_type = ["grf", SD, r_t_inactv_GSK3b])
    
    pn.add_transition(transition_id = 't_phos_tau',
                    label      =     "Phosphorylation of tau",
                    input_place_ids         =  ['p_tau', 'p_GSK3b_act', 'p_cas3'],
                    input_arc_weights  =  [1, 0, 0],
                    output_place_ids         =  ['p_tauP'],
                    output_arc_weights =  [1],
                    distribution_type = ["grf", SD, r_t_phos_tau])
    
    pn.add_transition(transition_id = 't_dephos_tau',
                    label      =     "Dephosphorylation of tau",
                    input_place_ids         =  ['p_tauP', 'p_Ca_cyto'],
                    input_arc_weights  =  [1, 0],
                    output_place_ids         =  ['p_tau'],
                    output_arc_weights =  [1],
                    distribution_type = ["grf", SD, r_t_dephos_tau])
            
    
    #Cholesterol homeostasis transitions
    
    pn.add_transition(transition_id = 't_LDLR_endocyto',
                    label      =     "LDLR endocytosis",
                    input_place_ids         =  ['p_ApoEchol_extra', 'p_chol_ER'],
                    input_arc_weights  =  [0, 0],
                    output_place_ids         =  ['p_chol_LE'],
                    output_arc_weights =  [354],
                    distribution_type = ["grf", SD, r_t_LDLR_endocyto]) 
    
    # pn.add_transition(transition_id = 't_chol_trans_LE_ER',
    #                 label      =     "Chol transport LE-ER",
    #                 input_place_ids         =  ['p_chol_LE'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_chol_ER'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_chol_trans_LE_ER])
    
    # pn.add_transition(transition_id = 't_chol_trans_LE_mito',
    #                 label      =     "Chol transport LE-mito",
    #                 input_place_ids         =  ['p_chol_LE'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_chol_mito'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_chol_trans_LE_mito])
    
    # pn.add_transition(transition_id = 't_chol_trans_LE_PM',
    #                 label      =     "Chol transport LE-PM",
    #                 input_place_ids         =  ['p_chol_LE'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_chol_PM'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_chol_trans_LE_PM])
    
    # pn.add_transition(transition_id = 't_chol_trans_PM_ER',
    #                 label      =     "Chol transport PM-ER",
    #                 input_place_ids         =  ['p_chol_PM'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_chol_ER'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_chol_trans_PM_ER])
    
    # pn.add_transition(transition_id = 't_chol_trans_ER_mito',
    #                 label      =     "Chol transport ER-mito",
    #                 input_place_ids         =  ['p_chol_ER'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_chol_mito'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_chol_trans_ER_mito])
    
    # pn.add_transition(transition_id = 't_CYP27A1_metab',
    #                 label      =     "Chol metab CYP27A1",
    #                 input_place_ids         =  ['p_chol_mito'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_27OHchol_intra'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_CYP27A1_metab])                                          
 
    # pn.add_transition(transition_id = 't_CYP11A1_metab',
    #                 label      =     "Chol metab CYP11A1",
    #                 input_place_ids         =  ['p_chol_mito'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_preg'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_CYP11A1_metab])  
                                                
    # pn.add_transition(transition_id = 't_CYP7B1_metab',
    #                 label      =     "27OHchol metab CYP7B1",
    #                 input_place_ids         =  ['p_27OHchol_intra'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_7HOCA'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_CYP7B1_metab])                                                  
    
    # pn.add_transition(transition_id = 't_27OHchol_endocyto',
    #                 label      =     "27OHchol endocyto",
    #                 input_place_ids         =  ['p_27OHchol_extra'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_27OHchol_extra', 'p_27OHchol_intra'],
    #                 output_arc_weights =  [1,1],
    #                 distribution_type = ["grf", 0.1, r_t_27OHchol_endocyto])  
                                                
    # pn.add_transition(transition_id = 't_CYP46A1_metab',
    #                 label      =     "Chol metab CYP46A1",
    #                 input_place_ids         =  ['p_chol_ER'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_24OHchol_intra'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_CYP46A1_metab])                                                  
     
    # pn.add_transition(transition_id = 't_24OHchol_exocyto',
    #                 label      =     "24OHchol exocyto",
    #                 input_place_ids         =  ['p_24OHchol_intra'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_24OHchol_extra'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0.1, r_t_24OHchol_exocyto])                                                  
      
    # pn.add_transition(transition_id = 't_chol_trans_PM_ECM',
    #                 label      =     "Chol transport PM-ECM",
    #                 input_place_ids         =  ['p_chol_PM', 'p_24OHchol_intra'],
    #                 input_arc_weights  =  [1, 0],
    #                 output_place_ids         =  [],
    #                 output_arc_weights =  [],
    #                 distribution_type = ["grf", 0.1, r_t_chol_trans_PM_ECM])                                                  
      
    # #ER retraction and collapse transitions 
    
    # pn.add_transition(transition_id = 't_RTN3_exp',
    #                     label = 'Expression rate of RTN3',
    #                     input_place_ids = [],
    #                     input_arc_weights = [],
    #                     output_place_ids = ['p_RTN3_PN'],
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_RTN3_exp])
    
    # pn.add_transition(transition_id = 't_LE_retro',
    #                     label = 'retrograde transport of LEs & ER',
    #                     input_place_ids = ['p_ATP','p_chol_LE','p_RTN3_axon', 'p_tau'],
    #                     input_arc_weights = [ATPcons_t_LE_trans, 0, 1, 0],                  
    #                     output_place_ids = ['p_ADP', 'p_RTN3_PN'],
    #                     output_arc_weights = [ATPcons_t_LE_trans, 1],
    #                     distribution_type = ["grf", 0.1, r_t_LE_retro])

    # pn.add_transition(transition_id = 't_LE_antero',
    #                     label = 'anterograde transport of LEs & ER',
    #                     input_place_ids = ['p_ATP','p_RTN3_PN', 'p_tau'], 
    #                     input_arc_weights = [ATPcons_t_LE_trans, 1, 0],
    #                     output_place_ids = ['p_ADP','p_RTN3_axon'],
    #                     output_arc_weights = [ATPcons_t_LE_trans, 1],
    #                     distribution_type = ["grf", 0.1, r_t_LE_antero]) 

    # pn.add_transition(transition_id = 't_RTN3_aggregation',
    #                     label = 'aggregation of monomeric RTN3 into HMW RTN3',
    #                     input_place_ids = ['p_RTN3_axon', 'p_RTN3_PN', 'p_Ab'], 
    #                     input_arc_weights = [1, 1, 0],
    #                     output_place_ids = ['p_RTN3_HMW_cyto'],
    #                     output_arc_weights = [1],
    #                     distribution_type=["grf", 0.1, r_t_RTN3_aggregation]) 

    # pn.add_transition(transition_id = 't_RTN3_auto',
    #                     label = 'functional autophagy of HMW RTN3',
    #                     input_place_ids = ['p_RTN3_HMW_cyto', 'p_RTN3_axon'], 
    #                     input_arc_weights = [1, 0],
    #                     output_place_ids = ['p_RTN3_HMW_auto'],
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_RTN3_auto]) 

    # pn.add_transition(transition_id = 't_RTN3_lyso',
    #                     label = 'functional delivery of HMW RTN3 to the lysosome',
    #                     input_place_ids = ['p_RTN3_HMW_auto', 'p_tau'], 
    #                     input_arc_weights = [1, 0],
    #                     output_place_ids = ['p_RTN3_HMW_lyso'],
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_RTN3_lyso]) 

    # pn.add_transition(transition_id = 't_RTN3_dys_auto',
    #                     label = 'dysfunctional autophagy of HMW RTN3',
    #                     input_place_ids = ['p_RTN3_HMW_cyto', 'p_RTN3_axon'], 
    #                     input_arc_weights = [1, 0],
    #                     output_place_ids = ['p_RTN3_HMW_dys1'],
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_RTN3_dys_auto])

    # pn.add_transition(transition_id = 't_RTN3_dys_lyso',
    #                     label = 'dysfunctional delivery of HMW RTN3 to the lysosome',
    #                     input_place_ids = ['p_RTN3_HMW_auto', 'p_RTN3_HMW_dys1', 'p_tau'], 
    #                     input_arc_weights = [1, 0, 0],
    #                     output_place_ids = ['p_RTN3_HMW_dys2'],
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_RTN3_dys_lyso])
    
    # # Energy metablism transitions
    # pn.add_transition(transition_id = 't_krebs', 
    #             label = 'Krebs cycle', 
    #             input_place_ids = ['p_ADP', 'p_Ca_mito', "p_Ab"],
    #             input_arc_weights = [1, 0, 0],
    #             output_place_ids = ['p_reduc_mito', 'p_ATP'], 
    #             output_arc_weights = [4, 1],
    #             distribution_type=["grf", 0.1, r_t_krebs])
    
    # pn.add_transition(transition_id = 't_ATP_hydro_mito', 
    #             label = 'ATP hydrolysis by cellular processes', 
    #             input_place_ids = ['p_ATP'],
    #             input_arc_weights = [1],
    #             output_place_ids = ['p_ADP'], 
    #             output_arc_weights = [1],
    #             distribution_type = ["grf", 0.1, r_t_ATP_hydro_mito])   
    
    # pn.add_transition(transition_id = 't_ETC', 
    #             label = 'Electron transport chain', 
    #             input_place_ids = ['p_reduc_mito', 'p_ADP', 'p_Ca_mito', 'p_ROS_mito', 'p_chol_mito', "p_Ab"],
    #             input_arc_weights = [22/3.96, 440, 0, 0, 0, 0],
    #             output_place_ids = ['p_ATP', 'p_ROS_mito'], 
    #             output_arc_weights = [440, 0.06],
    #             distribution_type = ["grf", 0.1, r_t_ETC])
    
    # pn.add_transition(transition_id = 't_ROS_metab', 
    #             label = 'Neutralization of ROS', 
    #             input_place_ids = ['p_ROS_mito', 'p_chol_mito'],
    #             input_arc_weights = [1, 0],
    #             output_place_ids = ['p_H2O_mito'], 
    #             output_arc_weights = [1],
    #             distribution_type = ["grf", 0.1, r_t_ROS_metab])    

    # # Output transitions: Cas3 for apoptosis
    # pn.add_transition(transition_id = 't_mito_dysfunc',
    #                     label = 'Mitochondrial complex 1 dysfunction',
    #                     input_place_ids = ['p_ROS_mito','p_Ab'],
    #                     input_arc_weights = [1, 0], 
    #                     output_place_ids = ['p_cas3'],         
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_mito_dysfunc])
    # # Cas3 inactivation
    # pn.add_transition(transition_id = 't_cas3_inact',
    #                     label = 'Caspase 3 inactivation',
    #                     input_place_ids = ['p_cas3'],
    #                     input_arc_weights = [1], 
    #                     output_place_ids = [],         
    #                     output_arc_weights = [],
    #                     distribution_type = ["grf", 0.1, r_t_cas3_inact])
    
    # pn.add_transition(transition_id = 't_ROS_gener_Ab',
    #                     label = 'ROS generation by Abeta',
    #                     input_place_ids = ['p_Ab'],
    #                     input_arc_weights = [0], 
    #                     output_place_ids = ["p_ROS_mito"],         
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_ROS_gener_Ab])
    
    
    
    # #Calcium homeostasis transitions 
    # pn.add_transition(transition_id = 't_Ca_imp',
    #                     label = 'VGCC/NMDA import channels',
    #                     input_place_ids = ['p_Ca_extra'],
    #                     input_arc_weights = [0],  
    #                     output_place_ids = ['p_Ca_cyto'],         
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_Ca_imp]) 

    # pn.add_transition(transition_id = 't_mCU',
    #                     label = 'Ca import into mitochondria via mCU',
    #                     input_place_ids = ['p_Ca_cyto', 'p_Ca_mito'],
    #                     input_arc_weights = [1,0], 
    #                     output_place_ids = ['p_Ca_mito'],         
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_mCU]) 

    # pn.add_transition(transition_id = 't_MAM',
    #                     label = 'Ca transport from ER to mitochondria',
    #                     input_place_ids = ['p_Ca_ER', 'p_Ca_mito'],
    #                     input_arc_weights = [1,0], 
    #                     output_place_ids = ['p_Ca_mito'],         
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_MAM]) 

    # pn.add_transition(transition_id = 't_RyR_IP3R',
    #                     label = 'Ca export from ER',
    #                     input_place_ids = ['p_Ca_extra', 'p_Ca_ER'],
    #                     input_arc_weights = [0,1], 
    #                     output_place_ids = ['p_Ca_cyto'],         
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_RyR_IP3R]) 

    # pn.add_transition(transition_id = 't_SERCA',
    #                     label = 'Ca import to ER',
    #                     input_place_ids = ['p_Ca_cyto','p_ATP'],
    #                     input_arc_weights = [1,0.5], 
    #                     output_place_ids = ['p_Ca_ER','p_ADP'],         
    #                     output_arc_weights = [1,0.5],
    #                     distribution_type = ["grf", 0.1, r_t_SERCA]) 

    # pn.add_transition(transition_id = 't_NCX_PMCA',
    #                     label = 'Ca efflux to extracellular space',
    #                     input_place_ids = ['p_Ca_cyto'],
    #                     input_arc_weights = [1],
    #                     output_place_ids = [],         
    #                     output_arc_weights = [],
    #                     switch_place_ids = ['p_on3'],
    #                     switch_arc_weights = [0],
    #                     distribution_type = ["grf", 0.1, r_t_NCX_PMCA])
    
    # pn.add_transition(transition_id = 't_mNCLX',
    #                     label = 'Ca export from mitochondria via mNCLX',
    #                     input_place_ids = ['p_Ca_mito'],
    #                     input_arc_weights = [1], 
    #                     output_place_ids = ['p_Ca_cyto'],         
    #                     output_arc_weights = [1],
    #                     distribution_type = ["grf", 0.1, r_t_mNCLX]) 

    # # Discrete on/of-switches calcium pacemaking
    # pn.add_transition(transition_id = 't_A',
    #                     label = 'A',
    #                     input_place_ids = ['p_on4'],
    #                     firing_condition = lambda a: a['p_on4']==1,
    #                     reaction_speed_function = lambda a: 1,
    #                     consumption_coefficients = [1], 
    #                     output_place_ids = ['p_Ca_extra'],         
    #                     production_coefficients = [1]) 
    
    # hfpn.add_transition_with_speed_function(
    #                     transition_id = 't_B',
    #                     label = 'B',
    #                     input_place_ids = ['p_Ca_extra'],
    #                     firing_condition = lambda a: a['p_Ca_extra']==1,
    #                     reaction_speed_function = lambda a: 1,
    #                     consumption_coefficients = [1], 
    #                     output_place_ids = ['p_on2'],         
    #                     production_coefficients = [1],
    #                     delay=0.5) 
    
    # hfpn.add_transition_with_speed_function(
    #                     transition_id = 't_C',
    #                     label = 'C',
    #                     input_place_ids = ['p_on2'],
    #                     firing_condition = lambda a: a['p_on2']==1,
    #                     reaction_speed_function = lambda a: 1,
    #                     consumption_coefficients = [1], 
    #                     output_place_ids = ['p_on3'],         
    #                     production_coefficients = [1],
    #                     delay=0) 

    # hfpn.add_transition_with_speed_function(
    #                     transition_id = 't_D',
    #                     label = 'D',
    #                     input_place_ids = ['p_on3'],
    #                     firing_condition = lambda a: a['p_on3']==1,
    #                     reaction_speed_function = lambda a: 1,
    #                     consumption_coefficients = [1], 
    #                     output_place_ids = ['p_on4'],         
    #                     production_coefficients = [1],
    #                     delay=0.5)

    # # Link to energy metabolism in that it needs ATP replenishment
    # hfpn.add_transition_with_mass_action(
    #                     transition_id = 't_NaK_ATPase',
    #                     label = 'NaK ATPase',
    #                     rate_constant =  k_t_NaK_ATPase,
    #                     input_place_ids = ['p_ATP', 'p_on3'],
    #                     firing_condition = lambda a: a['p_on3']==1,
    #                     consumption_coefficients = [1,0], 
    #                     output_place_ids = ['p_ADP'],         
    #                     production_coefficients = [1])
       
       # Run the network X times
    #a = {place.place_id:place.tokens for place in petri_net_model.places.values()}
    pn.run(1000, print_stats=False)
    
    #BSL: A good looking curve is, 2000 run steps, standard deviation of fixed 1 token, 200 starting tokens for asec. 10% of the mean gives a very smooth curve.

    # Plot the time-evolution of the system
    #input the place ids into this list for plotting
    list_for_plot = ['p_chol_LE'] 
    
    pn.plot_time_evolution(list_for_plot)


    # Generate block diagram of the Petri net

    #pn.generate_diagram("/Users/brand/Documents/Github/PN_Alzheimers_Parkinsons/blockdiags/AD")

if __name__ == "__main__":
    main()
