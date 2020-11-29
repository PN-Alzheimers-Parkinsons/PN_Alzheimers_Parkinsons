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
sys.path.insert(0, cwd[:(cwd.index(root_folder)+len(root_folder))] + os.sep + "Stochasticity-Net" + os.sep + "PDutils" + os.sep)

from Stochastic_PN_Architecture_v2 import *
from Brandon_PD_SN_parameters import *
from Brandon_PD_SN_rate_functions import *
from Brandon_PD_SN_initial_tokens import *
from Brandon_PD_SN_firing_conditions import *
from Stochastic_Analysis import Analysis


def main():
        #only runs this chunk of code if running the file directly, not as an import

        # Initialize an empty Petri net
    pn = PetriNet(number_of_runs=1) # so when you assign PetriNet class to pn, you also assign PetriNetModel class to the object petri_net_model since they are linked that way. This occurs because pn goes via "self" and what actually happens, is that you get pn.petri_net_model as an attribute?(right word?s)
    
#         # Add places for each chemical species
 #0.1 means 10%, 1 means 100% Standard deviation
    #relatieve standard deviation (10%, distribution and standard) #add michaeles menten
   
   #  # Cholesterol homeostasis
    # pn.add_place(it_p_chol_PM, place_id="p_chol_PM",label="Chol - perinuclear region")
    pn.add_place(it_p_chol_LE, place_id="p_chol_LE",label="Chol - late endosome")
    pn.add_place(it_p_chol_ER, place_id="p_chol_ER",label="Chol - ER")
    # pn.add_place(it_p_chol_mito, place_id="p_chol_mito",label="Chol - mitochondria")
    # pn.add_place(it_p_27OHchol_extra, place_id="p_27OHchol_extra",label="27-OH chol - extracellular")
    # pn.add_place(it_p_27OHchol_intra, place_id="p_27OHchol_intra",label="27-OH chol - intracellular")
    pn.add_place(it_p_ApoEchol_extra, place_id="p_ApoEchol_extra",label="ApoE - extracellular")
    pn.add_place(it_p_ApoEchol_EE, place_id="p_ApoEchol_EE",label="ApoE - Early endosome")
    # pn.add_place(it_p_7HOCA, place_id="p_7HOCA",label="7-HOCA")
    # pn.add_place(it_p_preg, place_id="p_preg",label="Pregnenolon")
    # pn.add_place(it_p_24OHchol_extra, place_id="p_24OHchol_extra",label="24OHchol extra")
    # pn.add_place(it_p_24OHchol_intra, place_id="p_24OHchol_intra",label="24OHchol intra")
  
   #  # Late endosome pathology 
   #  pn.add_place(it_p_LRRK2_mut, 'p_LRRK2_mut','LRRK2 - mutated')
   #  # PD specific places in cholesterol homeostasis
   #  pn.add_place(it_p_24OHchol_intra, "p_GBA1","GBA1")
   #  pn.add_place(it_p_SNCA_act_extra, "p_SNCA_act_extra","a-synuclein - extracellular")
   #  pn.add_place(it_p_SNCAApoEchol_extra, "p_SNCAApoEchol_extra","a-synuclein-ApoE complex - extracellular")
   #  pn.add_place(it_p_SNCAApoEchol_intra, "p_SNCAApoEchol_intra","a-synuclein-ApoE complex - intracellular")    
    
   #  # Lewy bodies
   #  pn.add_place(it_p_SNCA_act, "p_SNCA_act","SNCA - active")
   #  pn.add_place(it_p_VPS35, "p_VPS35", "VPS35")
   #  pn.add_place(it_p_SNCA_inact, "p_SNCA_inact", "SNCA - inactive")
   #  pn.add_place(it_p_SNCA_olig, "p_SNCA_olig", "SNCA - Oligomerised")
    pn.add_place(it_p_LB, "p_LB", "Lewy body")
   #  pn.add_place(it_p_Fe2, "p_Fe2", "Fe2 iron pool")
   #  # Energy metabolism
    # pn.add_place(it_p_ROS_mito, "p_ROS_mito", "ROS - mitochondria")
    # pn.add_place(it_p_H2O_mito, "p_H2O_mito", "H2O - mitochondria")
   #  pn.add_place(it_p_reduc_mito, "p_reduc_mito", "Reducing agents - mitochondria")
    # pn.add_place(it_p_cas3, "p_cas3","caspase 3 - mitochondria")
    # pn.add_place(it_p_DJ1, "p_DJ1","DJ1 mutant")    
    
   # # Drug places 
   #  pn.add_place(it_p_NPT200, place_id="p_NPT200", label = "Drug NPT200")
   #  pn.add_place(it_p_DNL151, place_id="p_DNL151", label = "Drug DNL151")
   #  pn.add_place(it_p_LAMP2A, place_id="p_LAMP2A", label = "Drug LAMP2A")
    
   #  # Calcium homeostasis
   #  pn.add_place(it_p_Ca_cyto, "p_Ca_cyto", "Ca - cytosole")
   #  pn.add_place(it_p_Ca_mito, "p_Ca_mito","Ca - mitochondria")
   #  pn.add_place(it_p_Ca_ER, "p_Ca_ER", "Ca - ER")
   #  pn.add_place(it_p_ADP, "p_ADP","ADP - Calcium ER import")
   #  pn.add_place(it_p_ATP, "p_ATP","ATP - Calcium ER import")    
    
   #  # Discrete on/of-switches calcium pacemaking

    # pn.add_place(500, "p_Ca_extra", "on1 - Ca - extracellular")
    # pn.add_place(0, "p_on2","on2")
    # pn.add_place(0, "p_on3","on3")
    # pn.add_place(0, "p_on4","on4")    
    
   #  # HMW RTN3 (cycling between different cellular compartments)
   #  pn.add_place(it_p_RTN3_HMW_cyto, place_id="p_RTN3_HMW_cyto", label="HMW RTN3 (cytosol)")
   #  pn.add_place(it_p_RTN3_HMW_auto, place_id="p_RTN3_HMW_auto", label="HMW RTN3 (autophagosome)")
   #  pn.add_place(it_p_RTN3_HMW_lyso, place_id="p_RTN3_HMW_lyso", label="HMW RTN3 (degraded in lysosome)")
   #  pn.add_place(it_p_RTN3_HMW_dys1, place_id="p_RTN3_HMW_dys1", label="HMW RTN3 (type I/III dystrophic neurites)")
   #  pn.add_place(it_p_RTN3_HMW_dys2, place_id="p_RTN3_HMW_dys2", label="HMW RTN3 (type II dystrophic neurites)")    
        
   #  # Two places that are NOT part of this subpathway, but are temporarily added for establishing proper connections
   #  # They will be removed upon merging of subpathways
   #  pn.add_place(it_p_tau, place_id="p_tau", label = "Unphosphorylated tau")
   #  pn.add_place(it_p_tauP, place_id="p_tauP", label = "Phosphorylated tau")
    
   #  # Monomeric RTN3 (cycling between axonal and perinuclear regions)
   #  pn.add_place(it_p_RTN3_axon, place_id="p_RTN3_axon", label="monomeric RTN3 (axonal)")
   #  pn.add_place(it_p_RTN3_PN, place_id="p_RTN3_PN", label="monomeric RTN3 (perinuclear)")
    #Cholesterol Homeostasis Transitions
    
    # Cholesterol Endocytosis
    # pn.add_transition(transition_id = 't_LDLR_endocyto',
    #                 label      =     "LDLR endocyto",
    #                 input_place_ids         =  ['p_ApoEchol_extra', 'p_chol_ER', 'p_LB'],
    #                 input_arc_weights  =  [0,0,0],
    #                 output_place_ids         =  ['p_ApoEchol_EE'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf", 0, r_t_LDLR_endocyto,fc_t_LDLR_endocyto ])
    
#     # Cleavage of cholesteryl esters 
    pn.add_transition(transition_id = 't_ApoEchol_cleav',
                    label      =     "ApoE-chol cleav",
                    input_place_ids         =  ['p_ApoEchol_EE'],
                    input_arc_weights  =  [1],
                    output_place_ids         =  ['p_chol_LE'],
                    output_arc_weights =  [354], 
                    distribution_type = ["grf",0,r_t_ApoEchol_cleav, fc_t_ApoEchol_cleav]) #
# #     # Transport Cholesterol from LE to ER
    # pn.add_transition(transition_id = 't_chol_trans_LE_ER',
    #                 label      =     "Chol transport LE-ER",
    #                 input_place_ids         =  ['p_chol_LE'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_chol_ER'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf",0,r_t_chol_trans_LE_ER, fc_t_chol_trans_LE_ER])
    
#     # Transport Cholesterol from LE to mito
    # pn.add_transition(transition_id = 't_chol_trans_LE_mito',
    #                 label      =     "Chol transport LE-mito",
    #                 input_place_ids         =  ['p_chol_LE'],
    #                 input_arc_weights  =  [1],
    #                 output_place_ids         =  ['p_chol_mito'],
    #                 output_arc_weights =  [1],
    #                 distribution_type = ["grf",0,r_t_chol_trans_LE_mito, fc_t_chol_trans_LE_mito])

# #     #Transport Cholesterol from LE to PM
#     pn.add_transition(transition_id = 't_chol_trans_LE_PM',
#                     label      =     "Chol transport LE-PM",
#                     input_place_ids         =  ['p_chol_LE'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_chol_PM'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,r_t_chol_trans_LE_PM,fc_t_chol_trans_LE_PM])
    
    
# #     # Transport Cholesterol from PM to ER
#     pn.add_transition(transition_id = 't_chol_trans_PM_ER',
#                     label      =     "Chol transport PM-ER",
#                     input_place_ids         =  ['p_chol_PM'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_chol_ER'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,r_t_chol_trans_PM_ER,fc_t_chol_trans_PM_ER])    
    


# #     # Transport Cholesterol from ER to PM
#     pn.add_transition(transition_id = 't_chol_trans_ER_PM',
#                     label      =     "Chol transport ER-PM",
#                     input_place_ids         =  ['p_chol_ER'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_chol_PM'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0,r_t_chol_trans_ER_PM, fc_t_chol_trans_ER_PM])   
    


# #     # Transport Cholesterol from ER to mito
#     pn.add_transition(transition_id = 't_chol_trans_ER_mito',
#                     label      =     "Chol transport ER-mito",
#                     input_place_ids         =  ['p_chol_ER'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_chol_mito'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,r_t_chol_trans_ER_mito,fc_t_chol_trans_ER_mito])
    

# #     # Metabolisation of chol by CYP27A1
#     pn.add_transition(transition_id = 't_CYP27A1_metab',
#                     label      =     "Chol metab CYP27A1",
#                     input_place_ids         =  ['p_chol_mito'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_27OHchol_intra'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,proper_rate_t_CYP27A1_metab,fc_t_CYP27A1_metab])                                          
    
    
# #     # Metabolism of chol by CYP11A1
#     pn.add_transition(transition_id = 't_CYP11A1_metab',
#                     label      =     "Chol metab CYP11A1",
#                     input_place_ids         =  ['p_chol_mito'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_preg'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,proper_rate_t_CYP11A1_metab, fc_t_CYP11A1_metab])  

# #     # Metabolisation of 27OHchol by CYP7B1
#     pn.add_transition(transition_id = 't_CYP7B1_metab',
#                     label      =     "27OHchol metab CYP7B1",
#                     input_place_ids         =  ['p_27OHchol_intra'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_7HOCA'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,proper_rate_t_CYP7B1_metab,fc_t_CYP7B1_metab])                                                  
    
    


# #     # Endocytosis of 27OHchol
#     pn.add_transition(transition_id = 't_27OHchol_endocyto',
#                     label      =     "27OHchol endocyto",
#                     input_place_ids         =  ['p_27OHchol_extra'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_27OHchol_extra', 'p_27OHchol_intra'],
#                     output_arc_weights =  [1,1], 
#                     distribution_type = ["grf",0,r_t_27OHchol_endocyto,fc_t_27OHchol_endocyto])  
    


# #     # Metabolisation of chol by CYP46A1
#     pn.add_transition(transition_id = 't_CYP46A1_metab',
#                     label      =     "Chol metab CYP46A1",
#                     input_place_ids         =  ['p_chol_ER'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_24OHchol_intra'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,proper_rate_t_CYP46A1_metab,fc_t_CYP46A1_metab])                                                  
    

# #     # Exocytosis of 24OHchol
#     pn.add_transition(transition_id = 't_24OHchol_exocyto',
#                     label      =     "24OHchol exocyto",
#                     input_place_ids         =  ['p_24OHchol_intra'],
#                     input_arc_weights  =  [1],
#                     output_place_ids         =  ['p_24OHchol_extra'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0,r_t_24OHchol_exocyto,fc_t_24OHchol_exocyto])                                                  
    

# #     # Transport of Chol into ECM
#     pn.add_transition(transition_id = 't_chol_trans_PM_ECM',
#                     label      =     "Chol transport PM-ECM",
#                     input_place_ids         =  ['p_chol_PM', 'p_24OHchol_intra'],
#                     input_arc_weights  =  [1, 0],
#                     output_place_ids         =  [],
#                     output_arc_weights =  [], 
#                     distribution_type = ["grf",0,r_t_chol_trans_PM_ECM, fc_t_chol_trans_PM_ECM])                                                  
    

# #resume
#     # Lewy bodies pathology
#     pn.add_transition(transition_id = 't_SNCA_degr',
#                     label      =     "SNCA degradation by CMA",
#                     input_place_ids         =  ['p_SNCA_act','p_VPS35','p_LRRK2_mut','p_27OHchol_intra','p_DJ1', 'p_DNL151', 'p_LAMP2A'],
#                     input_arc_weights  =  [1,0,0,0,0,0,0],
#                     output_place_ids         =  ['p_SNCA_inact'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_SNCA_degr])                                                  
    

    
#     pn.add_transition(transition_id = 't_SNCA_aggr',
#                     label      =     "SNCA aggregation",
#                     input_place_ids         =  ['p_SNCA_act','p_Ca_cyto','p_ROS_mito', 'p_tauP', 'p_NPT200'],
#                     input_arc_weights  =  [30,0,0,0,0], #should be reviewed if Ca is consumed
#                     output_place_ids         =  ['p_SNCA_olig'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_SNCA_aggr])                                                  



#     pn.add_transition(transition_id = 't_SNCA_fibril',
#                     label      =     "SNCA fibrillation",
#                     input_place_ids         =  ['p_SNCA_olig'],
#                     input_arc_weights  =  [100],
#                     output_place_ids         =  ['p_LB'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_SNCA_fibril])                                                  


#     pn.add_transition(transition_id = 't_IRE',
#                     label      =     "IRE",
#                     input_place_ids         =  ['p_Fe2'],
#                     input_arc_weights  =  [0],
#                     output_place_ids         =  ['p_SNCA_act'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_IRE]) #I think this catalyses?

#     # # PD specific
    
#     pn.add_transition(transition_id = 't_SNCA_bind_ApoEchol_extra',
#                     label      =     "Extracellular binding of SNCA to chol",
#                     input_place_ids         =  ['p_ApoEchol_extra','p_SNCA_act'],
#                     input_arc_weights  =  [0,30],
#                     output_place_ids         =  ['p_SNCA_olig'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_SNCA_bind_ApoEchol_extra])
    
#     pn.add_transition(transition_id = 't_chol_LE_upreg',
#                     label      =     "Upregulation of chol in LE",
#                     input_place_ids         =  ['p_GBA1'],
#                     input_arc_weights  =  [0],
#                     output_place_ids         =  ['p_chol_LE'],
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_chol_LE_upreg]) #BSL: I think this should be a catalysis Arc
    
#     # # Calcium homeostasis
#     pn.add_transition(transition_id = 't_Ca_imp',
#                     label      =     "L-type Ca channel",
#                     input_place_ids         =  ['p_Ca_extra'],
#                     input_arc_weights  =  [0],  # Need to review this 
#                     output_place_ids         =  ['p_Ca_cyto'],
#                     output_arc_weights =  [1], # Need to review this 
#                     distribution_type = ["grf",0.1,r_t_Ca_imp])
                    
#     pn.add_transition(transition_id = 't_mCU',
#                     label      =     "Ca import into mitochondria via mCU",
#                     input_place_ids         =  ['p_Ca_cyto','p_Ca_mito'],
#                     input_arc_weights  =  [1,0],  
#                     output_place_ids         =  ['p_Ca_mito'], 
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_mCU])

#     pn.add_transition(transition_id = 't_MAM',
#                     label      =     "Ca transport from ER to mitochondria",
#                     input_place_ids         =  ['p_Ca_ER','p_Ca_mito'],
#                     input_arc_weights  =  [1,0],  
#                     output_place_ids         =  ['p_Ca_mito'],  
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_MAM])

#     pn.add_transition(transition_id = 't_RyR_IP3R',
#                     label      =     "Ca export from ER",
#                     input_place_ids         =  ['p_Ca_extra','p_Ca_ER'],
#                     input_arc_weights  =  [0,1],  
#                     output_place_ids         =  ['p_Ca_cyto'], 
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_RyR_IP3R])
    
#     pn.add_transition(transition_id = 't_SERCA',
#                     label      =     "Ca import to ER",
#                     input_place_ids         =  ['p_Ca_cyto','p_ATP'],
#                     input_arc_weights  =  [1,1],  #!!!!! Need to review this 0 should be 1
#                     output_place_ids         =  ['p_Ca_ER','p_ADP'], 
#                     output_arc_weights =  [1,1], # Need to review this
#                     distribution_type = ["grf",0.1,r_t_SERCA])

#     pn.add_transition(transition_id = 't_NCX_PMCA',
#                     label      =     "Ca efflux to extracellular space",
#                     input_place_ids         =  ['p_Ca_cyto','p_on3'],
#                     input_arc_weights  =  [1,0],
#                     output_place_ids         =  [], 
#                     output_arc_weights =  [], 
#                     distribution_type = ["grf",0.1,r_t_NCX_PMCA])
    
#     pn.add_transition(transition_id = 't_mNCLX',
#                     label      =     "Ca export from mitochondria via mNCLX",
#                     input_place_ids         =  ['p_Ca_mito','p_LRRK2_mut'],
#                     input_arc_weights  =  [1,0], 
#                     output_place_ids         = ['p_Ca_cyto'], 
#                     output_arc_weights =   [1], 
#                     distribution_type = ["grf",0.1,r_t_mNCLX])
    
#     # Discrete on/of-switches calcium pacemaking
    
    # pn.add_transition(transition_id = 't_A',
    #                 label      =     "A",
    #                 input_place_ids         =  ['p_on4'],
    #                 input_arc_weights  =  [500],  
    #                 output_place_ids         = ['p_Ca_extra'], 
    #                 output_arc_weights =  [500], 
    #                 distribution_type = ["calcium",0.1,rate_t_A, fc_t_A]) #in HFPN, there is a delay of 0.5 s, maybe we need to be careful of this?

    # pn.add_transition(transition_id = 't_B',
    #                 label      =     "B",
    #                 input_place_ids         =  ['p_Ca_extra'],
    #                 input_arc_weights  =  [1],  
    #                 output_place_ids         = ['p_on2'], 
    #                 output_arc_weights =  [1], 
    #                 distribution_type = ["calcium",0.1,rate_t_B, fc_t_B])#in HFPN, there is a delay of 0.5 s, maybe we need to be careful of this?
    
    # pn.add_transition(transition_id = 't_C',
    #                 label      =     "C",
    #                 input_place_ids         =  ['p_on2'],
    #                 input_arc_weights  =  [500],  
    #                 output_place_ids         = ['p_on3'], 
    #                 output_arc_weights =  [500], 
    #                 distribution_type = ["calcium",0.1,rate_t_C, fc_t_C])

    # pn.add_transition(transition_id = 't_D',
    #                 label      =     "D",
    #                 input_place_ids         =  ['p_on3'],
    #                 input_arc_weights  =  [1],  
    #                 output_place_ids         = ['p_on4'], 
    #                 output_arc_weights =  [1], 
    #                 distribution_type = ["calcium",0.1,rate_t_D, fc_t_D])


    
#     # # Link to energy metabolism in that it needs ATP replenishment
#     #massaction problem
#     pn.add_transition(transition_id = 't_NaK_ATPase',
#                     label      =     "NaK ATPase",
#                     input_place_ids         =  ['p_ATP', 'p_on3'],
#                     input_arc_weights  =  [1,0],   
#                     output_place_ids         = ['p_ADP'], 
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,rate_t_NAK_ATPase])    
    
#     # # Energy metabolism
#     pn.add_transition(transition_id = 't_ATP_hydro_mito',
#                     label      =     "ATP hydrolysis in mitochondria",
#                     input_place_ids         =  ['p_ATP'],
#                     input_arc_weights  =  [1],    
#                     output_place_ids         = ['p_ADP'], 
#                     output_arc_weights =  [1], 
#                     distribution_type = ["grf",0.1,r_t_ATP_hydro_mito])       
    
    
    # pn.add_transition(transition_id = 't_ROS_metab',
    #                 label      =     "ROS neutralisation",
    #                 input_place_ids         =  ['p_ROS_mito','p_chol_mito','p_LB','p_DJ1'],
    #                 input_arc_weights  =  [1,0,0,0],    
    #                 output_place_ids         = ['p_H2O_mito'], 
    #                 output_arc_weights =  [1], 
    #                 distribution_type = ["grf",0.1,r_t_ROS_metab])   


#     # # #Link of krebs to calcium homeostasis
    
#     pn.add_transition(transition_id = 't_krebs',
#                     label      =     "Krebs cycle",
#                     input_place_ids         =  ['p_ADP','p_Ca_mito'],
#                     input_arc_weights  =  [1,0], #need review    
#                     output_place_ids         = ['p_reduc_mito','p_ATP'],         
#                     output_arc_weights =  [4,1], 
#                     distribution_type = ["grf",0.1,r_t_krebs])  

#     # #Link of ETC to calcium and cholesterol
    
#     pn.add_transition(transition_id = 't_ETC',
#                     label      =     "Electron transport chain",
#                     input_place_ids         =  ['p_reduc_mito', 'p_ADP', 'p_Ca_mito', 'p_chol_mito','p_ROS_mito','p_LRRK2_mut'],
#                     input_arc_weights  =  [22/3,22,0,0,0,0], # Need to review this   
#                     output_place_ids         = ['p_ATP', 'p_ROS_mito'],         
#                     output_arc_weights =  [22,0.005],
#                     distribution_type = ["grf",0.1,r_t_ETC])  
#     # # # Output transitions: Cas3 for apoptosis

    # pn.add_transition(transition_id = 't_mito_dysfunc',
    #                 label      =     "Mitochondrial complex 1 dysfunction",
    #                 input_place_ids         =  ['p_ROS_mito'],
    #                 input_arc_weights  =  [1],     
    #                 output_place_ids         = ['p_cas3'], 
    #                 output_arc_weights =  [1], 
    #                 distribution_type = ["grf",1,r_t_mito_dysfunc]) 
    
    # pn.add_transition(transition_id = 't_cas3_inact',
    #                 label      =     "Caspase 3 degredation",
    #                 input_place_ids         =  ['p_cas3'],
    #                 input_arc_weights  =  [1],     
    #                 output_place_ids         = [], 
    #                 output_arc_weights =  [], 
    #                 distribution_type = ["grf",1,r_t_cas3_inact,fc_t_cas3_inact]) 

    
#     # # Late endosome pathology #problem
#     pn.add_transition(transition_id = 't_phos_tau',
#                     label      =     "Phosphorylation of tau",
#                     input_place_ids         =  ['p_tau', 'p_SNCA_act'],
#                     input_arc_weights  =  [1, 0],     
#                     output_place_ids         = ['p_tauP'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,proper_rate_t_phos_tau])     
# #problem
#     pn.add_transition(transition_id = 't_dephos_tauP',
#                     label      =     "Dephosphorylation of tau protein",
#                     input_place_ids         =  ['p_tauP', 'p_Ca_cyto'],
#                     input_arc_weights  =  [1, 0],     
#                     output_place_ids         = ['p_tau'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,proper_rate_t_dephos_tauP]) 
    
#     pn.add_transition(transition_id = 't_RTN3_exp',
#                     label      =     "Expression rate of RTN3",
#                     input_place_ids         = [], 
#                     input_arc_weights  =  [],     
#                     output_place_ids         = ['p_RTN3_PN'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,r_t_RTN3_exp]) 

#     ATPcons_t_LE_trans = 0 #BSL: placeholder value coz idk
#     pn.add_transition(transition_id = 't_LE_retro',
#                     label      =     "retrograde transport of LEs & ER",
#                     input_place_ids         = ['p_ATP','p_chol_LE','p_RTN3_axon', 'p_tau','p_LRRK2_mut','p_LB'], 
#                     input_arc_weights  = [ATPcons_t_LE_trans, 0, 1, 0,0,0],     
#                     output_place_ids         = ['p_ADP','p_RTN3_PN'],
#                     output_arc_weights = [ATPcons_t_LE_trans, 1], 
#                     distribution_type = ["grf",0.1,r_t_LE_retro]) 
    
#     pn.add_transition(transition_id = 't_LE_antero',
#                     label      =     "anterograde transport of LEs & ER",
#                     input_place_ids         = ['p_ATP','p_RTN3_PN', 'p_tau'], # didn't connect p_tau yet
#                     input_arc_weights  =  [ATPcons_t_LE_trans, 1, 0], # tune these coefficients based on PD   
#                     output_place_ids         = ['p_ADP','p_RTN3_axon'],
#                     output_arc_weights =  [ATPcons_t_LE_trans, 1], # tune these coefficients based on PD
#                     distribution_type = ["grf",0.1,r_t_LE_antero]) # get later from NPCD

#     pn.add_transition(transition_id = 't_RTN3_aggregation',
#                     label      =     "aggregation of monomeric RTN3 into HMW RTN3",
#                     input_place_ids         = ['p_RTN3_axon', 'p_RTN3_PN'],  
#                     input_arc_weights  = [1, 1],     
#                     output_place_ids         = ['p_RTN3_HMW_cyto'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,r_t_RTN3_aggregation]) 
    
#     pn.add_transition(transition_id = 't_RTN3_auto',
#                     label      =     "functional autophagy of HMW RTN3",
#                     input_place_ids         = ['p_RTN3_HMW_cyto', 'p_RTN3_axon'],  
#                     input_arc_weights  =  [1, 0],     
#                     output_place_ids         = ['p_RTN3_HMW_auto'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,r_t_RTN3_auto])   
    
#     pn.add_transition(transition_id = 't_RTN3_lyso',
#                     label      =     "functional delivery of HMW RTN3 to the lysosome",
#                     input_place_ids         = ['p_RTN3_HMW_auto', 'p_tau'],  
#                     input_arc_weights  =  [1, 0],    
#                     output_place_ids         = ['p_RTN3_HMW_lyso'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,r_t_RTN3_lyso])  
    
#     pn.add_transition(transition_id = 't_RTN3_dys_auto',
#                     label      =     "dysfunctional autophagy of HMW RTN3",
#                     input_place_ids         = ['p_RTN3_HMW_cyto', 'p_RTN3_axon'], 
#                     input_arc_weights  =  [1, 0],    
#                     output_place_ids         = ['p_RTN3_HMW_dys1'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,r_t_RTN3_dys_auto])# tune later when data are incorporated

#     pn.add_transition(transition_id = 't_RTN3_dys_lyso',
#                     label      =     "dysfunctional delivery of HMW RTN3 to the lysosome",
#                     input_place_ids         = ['p_RTN3_HMW_auto', 'p_RTN3_HMW_dys1', 'p_tau'], 
#                     input_arc_weights  =  [1, 0, 0],    
#                     output_place_ids         = ['p_RTN3_HMW_dys2'],
#                     output_arc_weights =  [1],
#                     distribution_type = ["grf",0.1,r_t_RTN3_dys_lyso])
       
       # Run the network X times
    #a = {place.place_id:place.tokens for place in petri_net_model.places.values()}
    pn.run(100, print_stats=False)
    

    # Plot the time-evolution of the system
    #input the place ids into this list for plotting
    list_for_plot = ["p_ApoEchol_EE"] 
    
    pn.plot_time_evolution(list_for_plot)
    pn.timeseries_mean_for_place("p_ApoEchol_EE")
    analysis = Analysis(pn)
    run_save_name = "ggjiaogjioa"
    Analysis.store_to_file(analysis, run_save_name)
    print('Network saved to : "' + run_save_name+'.pkl"')

    # Generate block diagram of the Petri net

    #pn.generate_diagram("/Users/brand/Documents/Github/PN_Alzheimers_Parkinsons/blockdiags/BrandonPD")

if __name__ == "__main__":
    main()
