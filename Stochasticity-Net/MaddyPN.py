import random
import copy
import numpy as np
import matplotlib.pyplot as plt
# from blockdiag import parser, builder, drawer
#from blockdiag import parser, builder, drawer
# Examples: http://blockdiag.com/en/blockdiag/examples.html#simple-diagram

MAX_LABEL_CHARACTERS = 50

def check_label_length(label, limit = MAX_LABEL_CHARACTERS):
    if len(label) > limit:
            raise ValueError(f'Label \"{label}\" exceeds {MAX_LABEL_CHARACTERS} characters.')



class Place:
    def __init__(self, initial_tokens, place_id, label):
        """Put a place in the Petri net.
        
        Args:
            initial_tokens: numer of token that the place is initialized with
        """
        self.place_id = place_id
        self.label = label                # Label of the place
        self.tokens = initial_tokens

    def __str__(self):
        return f"{self.label} has {self.tokens} tokens"
    
    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, Place):
            return self.id == other.id
        return False
        
        # p_1 == p_2 -> true iff p_1.id == p_2.id

class BaseArc:
    def __init__(self, place, arc_weight):
        """Base class for making an arc in the Petri net.

        Args: 
            place: the one place acting as source/target of the arc
            arc_weight: the number of tokens removed/added from/to the place.
            
        """
        self.place = place
        self.arc_weight = arc_weight
        self.coefficient_scalar = 1


class InArc(BaseArc):  
    """An arc going from a place P to a transition T."""
    def __str__(self):
        return f"InArc from {self.place}"

    def fire(self):
        """Remove tokens from the input place."""
        self.place.tokens -= self.arc_weight*self.coefficient_scalar

    def allowed_firing(self):
        """Check whether the input place has enough tokens for the transition to occur."""
        return self.place.tokens >= self.arc_weight*self.coefficient_scalar

class OutArc(BaseArc):
    """An arc going from a transition T to a place P."""
    def fire(self):
        """Add tokens to the output place."""
        self.place.tokens += self.arc_weight*self.coefficient_scalar

    def __str__(self):
        return f"OutArc to {self.place}"

class InhibArc(BaseArc):  
    """An inhibitory arc going from a place P to a transition T."""
    def __str__(self):
        return f"InArc from {self.place}"

    def allowed_firing(self):
        """Arc weight is the number of tokens for which the transition is inhibited."""
        return self.place.tokens < self.arc_weight*self.coefficient_scalar

class CatalArc(BaseArc):
    """An arc in which tokens are not used but need to be present"""
    def __str__(self):
        return f"CatalArc from {self.place}"
    
    def allowed_firing(self):
        """Only fires if tokens are present"""
        return self.place.tokens >= self.arc_weight*self.coefficient_scalar #if no coefficient_scalar will produce a bug, where you still have firing, brandon
    
    #brandon 
    
    
class Transition:
    def __init__(self, in_arcs, out_arcs, inhib_arcs, catal_arcs, transition_id, label, distribution_type):
        """Put a transition T in the Petri net.
        
            Args:
                in_arcs (iterable): collection of ingoing arcs, to the transition
                out_arcs (iterable): collection of outgoing arcs, to the transition
                inhib_arcs (iterable): collection of inhibitory arcs, to the transition
                catal_arcs: collection of catalysis arcs, to the transition 
        """
        self.in_arcs = set(in_arcs)
        self.out_arcs = set(out_arcs)
        self.inhib_arcs = set(inhib_arcs)
        self.catal_arcs = set(catal_arcs)
        self.transition_id = transition_id
        self.label = label
        self.distribution_type = distribution_type
    def stochastic_distribution(self, distribution_type):
        
        pass


    def __str__(self):
        return f"Transition {self.label}"
        
    def fire(self):
        """Fire!"""
        
        #poisson and bernoulli try, gene transcribed maybe bernoulli
        #have the function in here? [g,2,1] # 10% standard deviation try, also try 20% SD, how does it differ for HFPN
        if self.distribution_type[0] == "g":
            random_integer = random.gauss(self.distribution_type[1],self.distribution_type[2])
        if self.distribution_type[0] =="u":
            random_integer = random.randint(self.distribution_type[1],self.distribution_type[2])


        #print("coefficient_scalar is")
        for arc2 in self.out_arcs.union(self.in_arcs):
            arc2.coefficient_scalar = random_integer
        # Check if transition can occur by inspecting whether each place has enough tokens
        firing_allowed = all(in_arc.allowed_firing() for in_arc in self.in_arcs)
        firing_not_inhibited = all(inhib_arc.allowed_firing() for inhib_arc in self.inhib_arcs)
        firing_not_inhibited2 = all(catal_arc.allowed_firing() for catal_arc in self.catal_arcs)
        if firing_allowed and firing_not_inhibited and firing_not_inhibited2:

            # Do "firing" for all input and output arcs associated with a transition
            for arc in self.out_arcs.union(self.in_arcs): 
                arc.fire()
        return firing_allowed and firing_not_inhibited and firing_not_inhibited2# Return if fired

class PetriNet:
    """Creates many equal copies of the specified Petri net for each run. 
    Has instance properties reporting on the statistical behaviour across all runs.
        timeseries_mean: mean number of tokens for each timestep (dimension 0) for each place (dimension 1)
        timeseries_std: var of tokens for each timestep (dimension 0) for each place (dimension 1)

    These statistical properties cannot be accessed before run()

    The initial conditions of the net is sotred in petri_net. This 
    original version of the net is not used for run operation. 
    """
    def __init__(self, number_of_runs):
        """Initializes an empty PetriNet 

            Args:
                number_of_runs (int): number of copies of the PetriNet to be run
        """
        self.locked = False     # initialized in state where editing allowed
        self.number_of_runs = number_of_runs
        self.petri_net_model = PetriNetModel()
        self.petri_net_copies = [] # for different runs

    def add_place(self, initial_tokens, place_id, label):
        """Add a place to the Petri net.
            Args: 
                initial_tokens (int): initial number of tokens in this place
                label (str): description of the place
            Return:
                place: an instance of the class Place, used to create transitions
        """

        if self.locked:
            raise RuntimeError('Error: do not change PetriNet after it has run.')

        self.petri_net_model.add_place(initial_tokens, place_id, label)

    def add_transition(self, transition_id, label, input_place_ids, input_arc_weights, output_place_ids, output_arc_weights, inhib_place_ids=[], inhib_arc_weights=[], catal_place_ids =[], catal_arc_weights=[], distribution_type=["u", 0, 1]):
        """Add a transition to the Petri net.
            Args:
                input_places (list): collection of places with arcs going into a transition
                output_places (list): collection of places with arcs going out of a transition
                inhib_places (list): collection of places with arcs inhibiting a transition
                catal_places (list): collection of places with arcs inhibiting a transition
        """
        if self.locked:
            raise RuntimeError('Error: do not change PetriNet after it has run.')

        self.petri_net_model.add_transition(transition_id, label, input_place_ids, input_arc_weights, output_place_ids, output_arc_weights, inhib_place_ids, inhib_arc_weights, catal_place_ids, catal_arc_weights, distribution_type)

    def run(self, number_of_steps, print_stats = False):
        """Runs multiple copies of the Petri net declared using this instance. 
        
        Args:
            number_of_steps (int): number of time-steps that the simulation will run for
            print_stats (bool): whether to print the mean and std of tokens across all copies/runs (for each time-step)
        """

        self.locked = True

        # Make deep copies of petri_net 
        self.petri_net_copies = [ self.petri_net_model.make_copy_of(self.petri_net_model) for _ in  range(self.number_of_runs) ]

        sorted_keys = self.petri_net_model.places.keys()
        self.place_id_to_index = dict(zip(sorted_keys, range(len(sorted_keys))))
        
        self.timeseries_mean = np.zeros((number_of_steps, len(self.place_id_to_index)))
        self.timeseries_std = np.zeros_like(self.timeseries_mean)
        
        # Run all copies and collect result (for later graphing for example)
        for step in range(number_of_steps):
            # runstep_tokens: array of arrays containing the number of tokens for a given run-step
            #     -> dimension 0: petri-net copy (in order of self.petri_net_copies)
            #     -> dimension 1: place (in order self.petri_net_model.places)
            
            runstep_tokens = [pn.run_one_step() for pn in self.petri_net_copies] 
            self.timeseries_mean[step] = np.mean(runstep_tokens, axis = 0) # averaging across the different copies for each place
            self.timeseries_std[step] = np.std(runstep_tokens, axis = 0) 

        if print_stats:
            print('Order of placs:\n', )
            print('Mean number of tokens for given time-step:\n', self.timeseries_mean)
            print('Var number of tokens for given time-step:\n', self.timeseries_std)

    def timeseries_mean_for_place(self, place_id):
        index = self.place_id_to_index[place_id]
        return self.timeseries_mean.T[index]

    def timeseries_std_for_place(self, place_id):
        index = self.place_id_to_index[place_id]
        return self.timeseries_std.T[index]


    def plot_time_evolution(self, place_ids = []):
        """Plots time-evolution using mean number of tokens for each place in the net.
        TODO: Only plot evolution of places spesified in place_ids if place_ids is not empty.
        """
#Info On the New Plotting Code Below:
        #places is a dictionary. So you can get the actual places using places.keys()
        #In the zip method below, self.petri_net_model.places.values() was changed to self.petri_net_model.places.keys() because this was more readable.
        
        dict_of_tokens = {} 
        for tokens, place in zip(self.timeseries_mean.T, self.petri_net_model.places.keys()):
            dict_of_tokens["{}".format(place)]=tokens #This line assigns each places list of tokens over all timesteps, to the dictionary dict_of_tokens and assigns the dictionary key using whats inside the square brackets.
        
        for x in place_ids:
            plt.plot(dict_of_tokens.get(x), label=x)        
        
        plt.legend(fontsize=20)
        plt.xlabel('Time-step')
        plt.ylabel('Mean tokens')
        plt.show()

    # def generate_diagram(self, file_name='diagram_petri_net'):
    #     petri_net_diagram = PetriNetDiagram(self.petri_net_model)
    #     petri_net_diagram.generate_image(file_name)


# Contains the information of a single petri net
# Designed for use by the PetriNet class, and should not be access from outside this file. 
class PetriNetModel:
    def __init__(self):
        """Initialize an empty Petri net."""
        self.places = {} # dictionary of [place_id : Place instance]
        self.transitions = {} # dictionary of [transition_id : Transition instance]
        self.successful_firings = []

    @staticmethod
    def make_copy_of(petri_net_model):
        """Makes a deep copy of a PetriNetModel instance. 
            
            Args: 
                petri_net_model: instance of PetriNetModel to be copied
        """

        pn_copy = PetriNetModel()
        for place in petri_net_model.places.values():
            pn_copy.add_place(place.tokens, place.place_id, place.label)

        # Copy over all transitions (referencing the newly created places)
        for t in petri_net_model.transitions.values():

            input_place_ids = [arc.place.place_id for arc in t.in_arcs]
            output_place_ids = [arc.place.place_id for arc in t.out_arcs]
            inhib_place_ids = [arc.place.place_id for arc in t.inhib_arcs]
            catal_place_ids = [arc.place.place_id for arc in t.catal_arcs]


            input_arc_weights = [arc.arc_weight for arc in t.in_arcs]
            output_arc_weights = [arc.arc_weight for arc in t.out_arcs]
            inhib_arc_weights = [arc.arc_weight for arc in t.inhib_arcs]
            catal_arc_weights = [arc.arc_weight for arc in t.catal_arcs]

            pn_copy.add_transition(    t.transition_id, t.label, input_place_ids, input_arc_weights, 
                                    output_place_ids, output_arc_weights, inhib_place_ids, inhib_arc_weights, catal_place_ids, catal_arc_weights, t.distribution_type)

        return pn_copy


    def add_place(self, initial_tokens, place_id, label):
        """Add a place to the Petri net. Prints warning if a place with the same id already exists. 
            
            Args: 
                initial_tokens (int): initial number of tokens in this place
                place_id (str): id of place
                label (str): description of the place
        """

        # Check if input is appropriate
        check_label_length(label)
        if ' ' in place_id:
            raise ValueError(f"Place_id should not contain any spaces {place_id}. Did you reverse place_id and label?")

        if place_id in self.places.keys():
            print(f"Warning: Place {place_id} already exists.")
        else:
            place = Place(initial_tokens, place_id, label)
            self.places[place_id] = place


    def add_transition(self, transition_id, label, input_place_ids, input_arc_weights, output_place_ids, output_arc_weights, inhib_place_ids=[], inhib_arc_weights=[], catal_place_ids=[], catal_arc_weights=[], distribution_type=[]):
        """Add a transition to the Petri net. Prints warning if a transition with the same id already exists. 
            Args:
                transition_id (str):
                label (str): 
                input_places (list): collection of places with arcs going into a transition
                output_places (list): collection of places with arcs going out of a transition
                inhib_places (list): collection of places with arcs inhibiting a transition
                catal_places (list): collection of places with arcs catalysing a transition
        """

        if len(input_place_ids) != len(input_arc_weights):
            raise ValueError(f"Unequal numbers of input places and input arc weights in transition {label}")
        if len(output_place_ids) != len(output_arc_weights):
            raise ValueError(f"Unequal numbers of output places and output arc weights in transition {label}")
        if len(inhib_place_ids) != len(inhib_arc_weights):
            raise ValueError(f"Unequal numbers of inhib places and inhib arc weights in transition {label}")
        if len(catal_place_ids) != len(catal_arc_weights):
            raise ValueError(f"Unequal numbers of cattal places and catal arc weights in transition {label}")

        check_label_length(label)
        if ' ' in transition_id:
            raise ValueError(f"Transition_id should not contain any spaces {transition_id}. Did you reverse transition_id and label?")

        if transition_id in self.transitions.keys():
            print(f"Warning: Transition {transition_id} already exists.")
        else:

            # translate input_places from strings to Place instances
            input_places = [self.places[place_id] for place_id in input_place_ids]
            output_places = [self.places[place_id] for place_id in output_place_ids]
            inhib_places = [self.places[place_id] for place_id in inhib_place_ids]
            catal_places = [self.places[place_id] for place_id in catal_place_ids]

            transition = Transition([InArc(place, weight) for (place, weight) in zip(input_places, input_arc_weights)],
                                    [OutArc(place, weight) for (place, weight) in zip(output_places, output_arc_weights)],
                                    [InhibArc(place, weight) for (place, weight) in zip(inhib_places, inhib_arc_weights)],
                                    [CatalArc(place, weight) for (place, weight) in zip(catal_places, catal_arc_weights)],
                                    transition_id, label, distribution_type)
            self.transitions[transition_id] = transition


    def run_one_step(self):
        """Run the Petri net for one step. 

            Returns: 
                A list of the current number of tokens for each place.
        """
        t = random.choice(list(self.transitions.values()))
        successful_firing = t.fire()
        self.successful_firings.append(successful_firing)

        return [place.tokens for place in self.places.values()]
    
# class PetriNetDiagram:

#       def __init__(self, petri_net):
#           """Initializes a Petri net diagram generator instance.

#             Args:
#                   petri_net: class instance of PetriNet to be plotted 
#           """
#           self.petri_net = petri_net


#       def build_string(self):
#           """Generates string needed for drawing diagram with blockdiag"""
#           strings = []

#           # List all places
#           for place_id, place in self.petri_net.places.items():
#               strings.append(f'{place_id} [label = "{place_id}\n{place.tokens}", shape = circle, color=lime];')

#           # List all transitions
#           for transition_id, transition in self.petri_net.transitions.items():
#               strings.append(f'{transition_id} [label = "{transition.transition_id}"];')

#               # Connect places to transitions (InArcs)
#               for in_arc in transition.in_arcs:
#                   strings.append(f'{in_arc.place.place_id} -> {transition_id} [label = "{in_arc.arc_weight}", fontsize = 11, color=black];')

#               # Connect transitions to places (OutArcs)
#               for out_arc in transition.out_arcs:
#                   strings.append(f'{transition_id} -> {out_arc.place.place_id} [label = "{out_arc.arc_weight}", fontsize = 11, color=blue];')
#             #connect catal arc brandon. I want to connect place to transition
#               for catal_arc in transition.catal_arcs:
#                   strings.append(f'{transition_id} -> {catal_arc.place.place_id} [label = "{catal_arc.arc_weight}", fontsize = 11, color=pink];')
                 
            
                 
                 

#           # build into string
#           return_string = 'blockdiag {\n'
#           for string in strings:
#               return_string += string + '\n'
#           return_string += '}'

#           return return_string

#           ## Example format of generated string
#           # blockdiag {
#           #    // Set labels to nodes.
#           #    A [label = "foo"];
#           #    B [label = "bar"];
#           #    // And set text-color
#           #    C [label = "baz"];

#           #    // Set labels to edges. (short text only)
#           #    A -> B [label = "click bar", textcolor="red"];
#           #    B -> C [label = "click baz"];
#           #    C -> A;
#           # }



#       def generate_image(self, file_name = 'petri-net-block-diagram'):
#           """Generates image file of block diagram.

#               Args:
#                   file_name (str): name of image-file
#           """

#           tree = parser.parse_string(self.build_string())
#           diagram = builder.ScreenNodeBuilder.build(tree)
#           draw = drawer.DiagramDraw('PNG', diagram, filename=file_name + '.png')
#           draw.draw()
#           draw.save()



if __name__ == "__main__":
        #only runs this chunk of code if running the file directly, not as an import

    
    #Initialize an empty Petri net
    pn = PetriNet(number_of_runs=100)
    
    #AB pathology
    pn.add_place(initial_tokens=1, place_id="p_asec", label="alpha secretase")
    pn.add_place(1, place_id="p_APP_PM", label="APP at plasma membrane")
    pn.add_place(0, place_id="p_APP_endo", label="endocytosed APP")
    pn.add_place(0, place_id="p_sAPPa", label="soluble sAPP alpha")
    pn.add_place(0, place_id="p_CTF83", label="CTF83")
    pn.add_place(1, place_id="p_bsec", label="beta secretase")
    pn.add_place(0, place_id="p_sAPPb", label="soluble sAPP beta")
    pn.add_place(0, place_id="p_CTF99", label="CTF99")
    pn.add_place(0, place_id="p_AB", label="Amyloid beta peptide")
    pn.add_place(0, place_id="p_AICD", label="AICD")
    pn.add_place(1, place_id="p_gsec", label="gamma secretase")

   
    #AB pathology transitions
    pn.add_transition(transition_id = 't_asec_exp',
                   label                = "alpha secretase expression",
                   input_place_ids       = [],
                   input_arc_weights  = [], 
                   output_place_ids       = ['p_asec'],
                   output_arc_weights = [1])
   
    pn.add_transition(transition_id = 't_asec_deg',
                   label      =     "alpha secretase degradation",
                   input_place_ids         =  ['p_asec'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  [],
                   output_arc_weights =  [])
   
    pn.add_transition(transition_id = 't_APP_asec_cleav',
                   label      =     "APP cleavage by alpha secretase",
                   input_place_ids         =  ['p_APP_PM'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_sAPPa', 'p_CTF83'],
                   output_arc_weights =  [1, 1],
                   catal_place_ids = ["p_asec"],
                   catal_arc_weights = [1]) #changed catal_arc_weight to 1 so that there is catalysis occurring brandon #also changing the catal_arc_weight doesnt change very much, oh wait, might need to add a new argument called catal_arc_weight, which will change the threshold of catalysis. Right now, the threshold of catalysis is determined by the ARC WEIGHT instead of the catal ARC WEIGHT
           
    pn.add_transition(transition_id = 't_APP_endo',
                   label      =     "APP endocytosis",
                   input_place_ids         =  ['p_APP_PM'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_APP_endo'],
                   output_arc_weights =  [1])
               
    pn.add_transition(transition_id = 't_APP_endo_bsec_cleav',
                   label      =     "APP cleavage by beta secretase",
                   input_place_ids         =  ['p_APP_endo'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_sAPPb', 'p_CTF99'],
                   output_arc_weights =  [1, 1],
                   catal_place_ids = ['p_bsec'],
                   catal_arc_weights= [1])
   
   
    pn.add_transition(transition_id = 't_bsec_exp',
                   label      =     "beta secretase expression",
                   input_place_ids         =  [],
                   input_arc_weights  =  [],
                   output_place_ids         =  ['p_bsec'],
                   output_arc_weights =  [1])
       
    pn.add_transition(transition_id = 't_bsec_deg',
                   label      =     "beta secretase degradation",
                   input_place_ids         =  ['p_bsec'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  [],
                   output_arc_weights =  [])
       
    pn.add_transition(transition_id = 't_CTF99_gsec_cleav',
                   label      =     "CTF99 cleavage by gamma secretase",
                   input_place_ids         =  ['p_CTF99'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_AB', 'p_AICD'],
                   output_arc_weights =  [1, 1],
                   catal_place_ids = ["p_gsec"],
                   catal_arc_weights = [1])
       
    pn.add_transition(transition_id = 't_gsec_exp',
                   label      =     "gamma secretase expression",
                   input_place_ids         =  [],
                   input_arc_weights  =  [],
                   output_place_ids         =  ['p_gsec'],
                   output_arc_weights =  [1])
       
    pn.add_transition(transition_id = 't_gsec_deg',
                   label      =     "gamma secretase degradation",
                   input_place_ids         =  ['p_gsec'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  [],
                   output_arc_weights =  [])
       
    pn.add_transition(transition_id = 't_bsec_deg',
                   label      =     "beta secretase degradation",
                   input_place_ids         =  ['p_bsec'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  [],
                   output_arc_weights =  [])
    
    
    #Tau pathology
    
    pn.add_place(initial_tokens=0, place_id="p_GSK3b_inact", label="Inactive GSK3 beta kinase")
    pn.add_place(1, place_id="p_GSK3b_act", label="Active GSK3 beta kinase")
    pn.add_place(0, place_id="p_tauP", label="Phosphorylated tau")
    pn.add_place(1, place_id="p_tau", label="Unphosphorylated tau (microtubule)")

    #Tau pathology transitions
    
    pn.add_transition(transition_id = 't_GSK3b_exp',
                   label                = "GSK3 beta expression",
                   input_place_ids       = [],
                   input_arc_weights  = [], 
                   output_place_ids       = ['p_GSK3b_inact'],
                   output_arc_weights = [1])
   
    pn.add_transition(transition_id = 't_GSK3b_deg',
                   label      =     "GSK3 beta degradation",
                   input_place_ids         =  ['p_GSK3b_inact'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  [],
                   output_arc_weights =  [])
   
    pn.add_transition(transition_id = 't_actv_GSK3b',
                   label      =     "GSK3 beta activation",
                   input_place_ids         =  ['p_GSK3b_inact'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_GSK3b_act'],
                   output_arc_weights =  [1])
                 
    pn.add_transition(transition_id = 't_inactv_GSK3b',
                   label      =     "GSK3 beta inactivation",
                   input_place_ids         =  ['p_GSK3b_act'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_GSK3b_inact'],
                   output_arc_weights =  [1])
    
    pn.add_transition(transition_id = 't_phos_tau',
                   label      =     "Phosphorylation of tau",
                   input_place_ids         =  ['p_tau'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_tauP'],
                   output_arc_weights =  [1],
                   catal_place_ids= ['p_GSK3b_act'],
                   catal_arc_weights= [1])
    
    pn.add_transition(transition_id = 't_dephos_tau',
                   label      =     "Dephosphorylation of tau",
                   input_place_ids         =  ['p_tauP'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_tau'],
                   output_arc_weights =  [1])
    
    pn.add_transition(transition_id = 't_inactv_GSK3b',
                   label      =     "GSK3 beta inactivation",
                   input_place_ids         =  ['p_GSK3b_act'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_GSK3b_inact'],
                   output_arc_weights =  [1])
            
    
    #Cholesterol homeostasis 
    # Cholesterol-ApoE
    pn.add_place(10, place_id="p_ApoEchol_extra", label="ApoE-cholesterol complex (extracellular)")
    #Cholesterol in different organelles
    pn.add_place(10, place_id="p_chol_LE", label="Cholesterol (late endosome)")
    pn.add_place(1, place_id="p_chol_mito", label="Cholesterol (mitochondria)")
    pn.add_place(1, place_id="p_chol_ER", label="Cholesterol (ER)")
    pn.add_place(1, place_id="p_chol_PM", label="Cholesterol (plasma membrane)") 
     # Oxysterols
    pn.add_place(1, place_id="p_24OHchol_extra", label="24-hydroxycholesterol (extracellular)")
    pn.add_place(1, place_id="p_24OHchol_intra", label="24-hydroxycholesterol (intracellular)")
    pn.add_place(1, place_id="p_27OHchol_extra", label="27-hydroxycholesterol (extracellular)")
    pn.add_place(1, place_id="p_27OHchol_intra", label="27-hydroxycholesterol (intracellular)")
    pn.add_place(1, place_id="p_7HOCA", label="7-HOCA")
    pn.add_place(1, place_id="p_preg", label="Pregnenolon")
    
    #Transitions
    pn.add_transition(transition_id = 't_LDLR_endocyto',
                   label      =     "LDLR endocytosis",
                   input_place_ids         =  ['p_ApoEchol_extra'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_LE'],
                   output_arc_weights =  [1],
                   inhib_place_ids= ['p_chol_ER'],
                   inhib_arc_weights= [1],
                   distribution_type=["g",5,1]) #not sure what to use for this 
    
    pn.add_transition(transition_id = 't_chol_trans_LE_ER',
                   label      =     "Chol transport LE-ER",
                   input_place_ids         =  ['p_chol_LE'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_ER'],
                   output_arc_weights =  [1])
    
    pn.add_transition(transition_id = 't_chol_trans_LE_mito',
                   label      =     "Chol transport LE-mito",
                   input_place_ids         =  ['p_chol_LE'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_mito'],
                   output_arc_weights =  [1])
    
    pn.add_transition(transition_id = 't_chol_trans_LE_PM',
                   label      =     "Chol transport LE-PM",
                   input_place_ids         =  ['p_chol_LE'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_PM'],
                   output_arc_weights =  [1])
    
    pn.add_transition(transition_id = 't_chol_trans_PM_ER',
                   label      =     "Chol transport PM-ER",
                   input_place_ids         =  ['p_chol_PM'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_ER'],
                   output_arc_weights =  [1])
    
    pn.add_transition(transition_id = 't_chol_trans_ER_mito',
                   label      =     "Chol transport ER-mito",
                   input_place_ids         =  ['p_chol_ER'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_chol_mito'],
                   output_arc_weights =  [1])
    
    pn.add_transition(transition_id = 't_CYP27A1_metab',
                   label      =     "Chol metab CYP27A1",
                   input_place_ids         =  ['p_chol_mito'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_27OHchol_intra'],
                   output_arc_weights =  [1])                                          
 
    pn.add_transition(transition_id = 't_CYP11A1_metab',
                   label      =     "Chol metab CYP11A1",
                   input_place_ids         =  ['p_chol_mito'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_preg'],
                   output_arc_weights =  [1])  
                                                
    pn.add_transition(transition_id = 't_CYP7B1_metab',
                   label      =     "27OHchol metab CYP7B1",
                   input_place_ids         =  ['p_27OHchol_intra'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_7HOCA'],
                   output_arc_weights =  [1])                                                  
    
    pn.add_transition(transition_id = 't_27OHchol_endocyto',
                   label      =     "27OHchol endocyto",
                   input_place_ids         =  ['p_27OHchol_extra'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_27OHchol_extra', 'p_27OHchol_intra'],
                   output_arc_weights =  [1,1])
    
    # This is an artificial transition !!
                                                
    pn.add_transition(transition_id = 't_CYP46A1_metab',
                   label      =     "Chol metab CYP46A1",
                   input_place_ids         =  ['p_chol_ER'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_24OHchol_intra'],
                   output_arc_weights =  [1])                                                  
     
    pn.add_transition(transition_id = 't_24OHchol_exocyto',
                   label      =     "24OHchol exocyto",
                   input_place_ids         =  ['p_24OHchol_intra'],
                   input_arc_weights  =  [1],
                   output_place_ids         =  ['p_24OHchol_extra'],
                   output_arc_weights =  [1])                                                  
      
    pn.add_transition(transition_id = 't_chol_trans_PM_ECM',
                   label      =     "Chol transport PM-ECM",
                   input_place_ids         =  ['p_chol_PM', 'p_24OHchol_intra'],
                   input_arc_weights  =  [1, 0],
                   output_place_ids         =  [],
                   output_arc_weights =  [])                                               
     

 
       # Run the network X times
    pn.run(100, print_stats=False)

    # Plot the time-evolution of the system
    places_to_plot = ['p_asec', 'p_tau', 'p_tauP']
    pn.plot_time_evolution(places_to_plot)


    # Generate block diagram of the Petri net
    # pn.generate_diagram("maddy_diagram")
