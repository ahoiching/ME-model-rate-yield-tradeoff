import cobrame
import json
from cobrame.io.json import load_reduced_json_me_model, load_json_me_model
#import pickle

with open('./iJL1678b.json', 'rb') as f:
    me = load_json_me_model(f)

#with open('./iJL1678b.pickle', 'rb') as f:
#    me = pickle.load(f)
    
def solve_me_model(me, max_mu, precision=1e-6, min_mu=0, using_soplex=True,
                  compiled_expressions=None):
    if using_soplex:
        from cobrame.solve.algorithms import binary_search
        binary_search(me, min_mu=min_mu, max_mu=max_mu, debug=True, mu_accuracy=precision,
                      compiled_expressions=compiled_expressions)
    else:
        from qminospy.me1 import ME_NLP1
        # The object containing solveME methods--composite that uses a ME model object 
        me_nlp = ME_NLP1(me, growth_key='mu')
        # Use bisection for now (until the NLP formulation is worked out)
        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=max_mu)
        me.solution.f = me.solution.x_dict['biomass_dilution']
        
solve_me_model(me, 1.5, min_mu = .1, precision=1e-6)

with open('s1_solution.json', 'w') as fp:
    json.dump(me.solution.x_dict, fp)

#print me.reactions.dummy_reaction_FWD_SPONT.objective_coefficient
