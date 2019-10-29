import cobrame
import json
from cobrame.io.json import load_reduced_json_me_model, load_json_me_model
from cobrame.solve.algorithms import binary_search,solve_at_growth_rate
import numpy as np
#import pickle

with open('./iJL1678b.json', 'rb') as f:
    me = load_json_me_model(f)

#with open('./iJL1678b.pickle', 'rb') as f:
#    me = pickle.load(f)

me.unmodeled_protein_fraction=0.18  #Changing the unmodeled fraction
me.gam=34.98
me.ngam=1.0

#TCA modification
TCA_reaction_id=['ACONTa','ACONTb','AKGDH','CS','FUM','ICDHyr','MDH','SUCDi','SUCOAS']
TCA_keff_dict={}
for r in me.reactions:
    if hasattr(r, 'stoichiometric_data'):
       #print(r.subsystem, ' ', r.stoichiometric_data.id)
       if r.stoichiometric_data.id in TCA_reaction_id and r.stoichiometric_data.id not in TCA_keff_dict:
           temp_dict=[]
           for r2 in r.stoichiometric_data.parent_reactions:
               temp_dict.append(r2.keff)
           TCA_keff_dict[r.stoichiometric_data.id]=temp_dict
print(TCA_keff_dict)

TCA_keff_dict_max={}
for tkd in TCA_keff_dict:
    TCA_keff_dict_max[tkd]=np.max(TCA_keff_dict[tkd])
    
MS_ATP_res=[37.2,30.22,21.46,15.22]
MS_data={
    'ACONTa':[0.92, 0.84, 0.66, 0.57], #acnA or acnB
    'ACONTb':[0.92, 0.84, 0.66, 0.57], #acnA or acnB
    'AKGDH':[0.33+0.38+0.63, 0.36+0.39+0.61, 0.29+0.35+0.56, 0.25+0.30+0.52], #sucA and sucB and Ipd
    'CS':[0.88, 0.80, 0.61, 0.48], #gltA
    'FUM':[0.24, 0.21, 0.17, 0.13], #fumB or fumC or fumA
    'ICDHyr':[1.55, 1.55, 1.31, 1.39],#icd
    'MDH':[0.45, 0.45, 0.41, 0.39],#mdh
    'SUCDi':[0.33+0.16+0+0, 0.3+0.15+0+0, 0.29+0.13+0+0, 0.24+0.11+0+0],#sdhA and sdhB and sdhC and sdhD
    'SUCOAS':[0.52+0.36, 0.49+0.35, 0.38+0.28, 0.30+0.22], #sucC and sucD
}

#linear regression:
x = np.array([0, 1, 2, 3])
y = np.array([-1, 0.2, 0.9, 2.1])

A = np.vstack([x, np.ones(len(x))]).T

m, c = np.linalg.lstsq(A, y)[0]

MS_data_mean={}
MS_data_linear_fit={}

y=np.array(MS_ATP_res)

for md in MS_data:
    MS_data_mean[md]=np.mean(MS_data[md])
    
    #linear regression
    x=np.array(MS_data[md])
    A = np.vstack([x, np.ones(len(x))]).T
    
    m, c = np.linalg.lstsq(A, y)[0]
    
    MS_data_linear_fit[md]=m
    
item=[]
value_model=[]
value_ms=[]

for mdm in MS_data_linear_fit:
    item.append(mdm)
    value_model.append(TCA_keff_dict_max[mdm])
    value_ms.append(MS_data_linear_fit[mdm])
    
#modifying the keffs
for mdlf in MS_data_linear_fit:
    for sd in me.stoichiometric_data:
        if sd.id in mdlf:
            for r in sd.parent_reactions:
                r.keff=r.keff*MS_data_linear_fit[mdlf]/TCA_keff_dict_max[mdlf]
                r.update()

TCA2_1=['AKGDH','SUCOAS','SUCDi']

div_coeff=4.0 #By iteration, it turns out that 3 lower stream TCA reactions need to tune down the protein effiencies even further.
for sd in me.stoichiometric_data:
    if sd.id in TCA2_1:
        for r in sd.parent_reactions:
            r.keff=r.keff/div_coeff
            r.update()
            
blk_met_rxn=[
#bp2 reactions:
#21
'EAR100x1', 'EAR40x1', 'EAR60x1', 'EAR80x1', 
#20    
'RNTR3c21',
 'EAR120x1',
 'RNTR4c21',
 'POR51',
 u'FACOAE100',
 'EAR140x1',
 u'FACOAE80',
 'RNTR1c21',
 'EAR141x1',
 'RNTR2c21',
 'EAR121x1',
 'AACPS81',
 'AACPS91',
 #19 
'AACPS21',
 'RNTR1c22',
 u'FACOAE140',
 u'FACOAE120',
 'AACPS71',
 u'CTECOAI6',
 'RNTR4c22',
 'RNTR3c22',
 'AACPS11',
 'RNTR2c22',
 'EAR160x1',
 u'FACOAE141',
 'EAR161x1',
 'EDD',#18
 #bp1 reactions
 #17
'AACPS41',
 u'R1PK',
 u'XPPT',
 u'HXAND',
 u'NTD11',
 u'FACOAE161',
 u'R15BPK',
 u'PUNP5',
 'AACPS31',
 u'CTECOAI7',
    
 'NADPHQR2',#16
 'VPAMTr',#15
 'PYAM5PO',#14
 'FLDR21',#13
 'FLDR22',#12
 'ALATA_L2',#11
 'TRSARr',#10
 'FTHFD',#9
 'DAAD',#8
 'FADRx2',#7
 'IDOND',#6
 'ASPT', #5
 'NADTRHD',#4
 'GLYAT', #3
 'ABTA', #2
 'ICL' #1
]

for sd in me.stoichiometric_data:
    if sd.id in blk_met_rxn:
        for r in sd.parent_reactions:
            r.upper_bound=0
            r.lower_bound=0
            
#blocking the exchange reactions
for r in me.reactions:
    if 'EX_' in r.id[0:3] and r.id not in ['EX_meoh_e' ,
'EX_so4_e' ,
'EX_RNase_m23' ,
'EX_RNase_m5' ,
'EX_nh4_e' ,
'EX_cu2_e' ,
'EX_co2_e' ,
'EX_pi_e' ,
'EX_5mtr_e' ,
'EX_mn2_e' ,
'EX_zn2_e' ,
'EX_slnt_e' ,
'EX_h_e' ,
'EX_cobalt2_e' ,
'EX_mg2_e' ,
'EX_mobd_e' ,
'EX_h2o_e' ,
'EX_tl_c' ,
'EX_fe2_e' ,
'EX_ac_e' ,
'EX_k_e' ,
'EX_o2_e' ,
'EX_glc__D_e' ,
'EX_RNase_m16']:
            r.upper_bound=0
            r.lower_bound=0
            
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
upf_list=[0.18,0.22,0.26,0.30]
final_dict={}
for upf in upf_list:
    print('Setting up unmodeled protein fraction...',upf)
    me.unmodeled_protein_fraction=upf  #Changing the unmodeled fraction
    compiled_expressions=compile_expressions(me)        
    solve_me_model(me, 1.5, min_mu = .1, precision=1e-2,compiled_expressions=compiled_expressions)

    #Get the maximum growth rate of the current model
    sol = me.solution
    muopt=me.solution.x_dict['biomass_dilution']

    #picking three different growth rates
    growth_list=[0.85*muopt,0.90*muopt,muopt]
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
            for direction in ['max']:
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
    final_dict[upf]=target_dict

with open('s4_solution_upf.json', 'w') as fp:
    json.dump(final_dict, fp)
