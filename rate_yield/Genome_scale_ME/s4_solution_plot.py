import numpy as np
import matplotlib.pyplot as plt
import json

with open("./s4_solution.json", 'r') as f:
    simulation_results = json.load(f)
growth_rates_str=np.array(list(simulation_results['SUR'].keys()))

ac_list=[]
glc_list_max=[]
glc_list_min=[]
growths=[]
for i in range(len(growth_rates_str)):
    growths.append(simulation_results['SUR'][growth_rates_str[i]]['max']['biomass_dilution'])
    ac_list.append(simulation_results['SUR'][growth_rates_str[i]]['max']['EX_ac_e']) #acetate excretion under maximum yield
    
    glc_list_max.append(simulation_results['SUR'][growth_rates_str[i]]['max']['EX_glc__D_e'])
    glc_list_min.append(simulation_results['SUR'][growth_rates_str[i]]['min']['EX_glc__D_e'])

ac_list=np.array(ac_list)
glc_list_max=np.array(glc_list_max)
glc_list_min=np.array(glc_list_min)
growths=np.array(growths)
    
#Experimental data
convert=0.5
MU_Basan = [0.88,0.81,0.75,0.71,0.68,0.58,0.97,0.94,0.91,0.88,0.87,0.81,0.78,0.74]

GUR_Basan = [5.58,4.94,4.31,3.9,3.94,3.25,5.4,4.97,5.48,5.11,5.15,4.69,4.46,4.21]
GUR_Basan = np.multiply(GUR_Basan, 1/convert)

AC_Basan = [1.88,0.7,0,0,0,0,2.52,2.37,2.38,1.94,1.5,1.16,0.03,0]
AC_Basan = np.multiply(AC_Basan, 1/convert)

GY_Basan = MU_Basan/GUR_Basan/180.1559 * 1000

plt.figure(figsize=(6,3))
plt.subplot(1,2,1)
plt.fill_between(growths,-growths/glc_list_max/180.1559 * 1000,-growths/glc_list_min/180.1559 * 1000,
                alpha=0.3)
plt.plot(growths,-growths/glc_list_max/180.1559 * 1000,'C0')
plt.plot(MU_Basan,GY_Basan,'C0^')
plt.ylabel('Growth yield(gDW/g)')
plt.xlabel('Growth(1/hr)')

plt.subplot(1,2,2)
plt.plot(growths,ac_list)
plt.plot(MU_Basan,AC_Basan,'C0^')
plt.ylabel('Excretion rate(mmol/gDW/hr)')
plt.xlabel('Growth(1/hr)')

plt.tight_layout()

plt.savefig('s4_solution_plot.pdf')
