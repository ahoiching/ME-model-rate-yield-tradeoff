import cobrame
import json
from cobrame.io.json import load_reduced_json_me_model, load_json_me_model
from cobrame.solve.algorithms import binary_search,solve_at_growth_rate
#import pickle

with open('./iJL1678b.json', 'rb') as f:
    me = load_json_me_model(f)

#with open('./iJL1678b.pickle', 'rb') as f:
#    me = pickle.load(f)
    
def solve_me_model(me, max_mu, precision=1e-2, min_mu=0, using_soplex=True,
                  compiled_expressions=None):
    if using_soplex:
        binary_search(me, min_mu=min_mu, max_mu=max_mu, debug=True, mu_accuracy=precision,
                      compiled_expressions=compiled_expressions)
    else:
        from qminospy.me1 import ME_NLP1
        # The object containing solveME methods--composite that uses a ME model object 
        me_nlp = ME_NLP1(me, growth_key='mu')
        # Use bisection for now (until the NLP formulation is worked out)
        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=max_mu)
        me.solution.f = me.solution.x_dict['biomass_dilution']

from cobrame.solve.symbolic import compile_expressions
compiled_expressions=compile_expressions(me)        
solve_me_model(me, 1.5, min_mu = .1, precision=1e-2,compiled_expressions=compiled_expressions)

#Get the maximum growth rate of the current model
sol = me.solution
muopt=me.solution.x_dict['biomass_dilution']

#picking three different growth rates
growth_list=[0.90*muopt, 0.95*muopt, muopt]
growth_list=sorted(growth_list,reverse=True)

#Now we start simulating the growth rate dependence 
#maximum and minimum glucose uptake (max and min growth yield)
target_dict={} #initializing the final solution dictionary
sur=me.reactions.EX_glc__D_e

if 'SUR' in target_dict:
    gr_sols_sur=target_dict['SUR']
elif 'SUR' not in target_dict:
    gr_sols_sur={}
for gr in sorted(growth_list):
    if gr<=muopt and gr not in gr_sols_sur:
        print('setting up gr=%.2f-----------------' %gr)
        sur.upper_bound=1000
        sur.lower_bound=-1000

        dir_sols={}   
        for direction in ['min','max']:
            me.objective={sur.id: 1.0 if direction is 'max' else -1.0}
            print('solving')
            solve_at_growth_rate(me,gr,compiled_expressions=compiled_expressions)
            status=me.solution.status
            print(gr,direction,me.solution.x_dict[sur.id],status)

            if status != 'optimal':
                dir_sols[direction]=None
            else:
                dir_sols[direction]=me.solution.x_dict
        gr_sols_sur[gr] = dir_sols
        target_dict['SUR']=gr_sols_sur

with open('s2_solution.json', 'w') as fp:
    json.dump(target_dict, fp)
