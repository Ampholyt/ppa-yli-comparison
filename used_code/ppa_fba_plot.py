# script should reporduce experimental data comparison with simulated data according to https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-11-83/tables/1
# imports
import cobra
import yaml
import numpy as np # plotting
import matplotlib.pyplot as plt # plotting
from cobra.flux_analysis import pfba
from cobra.medium import minimal_medium
import logging # remove the warnings while loading model => faster

import sys # append path

sys.path.append('../scripts/')
import helperFunction as hf

# config
config_name = 'model_config'
config_path = f'../config/{config_name}.yaml'

# load config
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

_verbose = False
_use_precomputed = False
_plotting = True
_without_iMT1026v1 = True # because iTM1026v1 has no meaningful results

## ppa: https://www.notion.so/Growth-rates-of-engineered-strains-fd69830a7b13471da65404a24c4a7abe?pvs=4#808d8f5f4f894c8e99baa3205455c664
# global variables (experimental data + maintenance)
experiment_number = 3 # set 3 to only simulate wild-type and ignore recombinant strain
glucose_uptake = [1.00, 1.28, 1.72, 1.01, 1.37, 1.56]
oxygen_uptake = [2.35, 2.01, 2.01, 2.44, 1.99, 1.81]
ppa_experimental_growth = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
ppa_co2_prod = [2.43, 2.55, 3.21, 2.52, 2.68, 2.94]
ppa_ethanol_prod = [1000, 0.31, 0.84, 1000, 0.41, 0.83]
# ppa_arabitol_prod = [xxxxx] # not used according to iMT1026v1 supplementary

# _ppa_maintenance = 2.26 # mmol/gDW/h # according to iLC915
_ppa_maintenance = 2.81 # ATP mmol/gDW/h # according to iMT1026

_condition_name = f'Maintenance set to {_ppa_maintenance} ATP mmol/gDW/h'

_growth_conditions = {'experimental': ppa_experimental_growth[:3]}
_co2_productions = {'experimental': ppa_co2_prod[:3]}
_condition_names = ['normoxic', 'O2-limited', 'hypoxic']



print('iLC915 model: ')
# load iLC915
iLC915_model = cobra.io.read_sbml_model(config['models']['ppaiLC915'])
# all relevant reactions:
iLC915_biomass_reaction = 'r1133' # biomass
iLC915_glu_ex_rxn = 'r1145' # glucose uptake
iLC915_fructose_ex_rxn = 'r1144' # fructose uptake
iLC915_glycerol_ex_rxn = 'r1148' # glycerol uptake
iLC915_o2_ex_rxn = 'r1160' # O2
iLC915_arabitol_ex_rxn = 'r1128' # Arabitol/arabinitol
iLC915_ethanol_ex_rxn = 'r1141' # ethanol
iLC915_maintenance = 'r1188' # maintenance 2.9 in iMT1026v3 (0 in iLC915)

# prepare iLC915: needed because of internal cycles and active reactions
# lower and upper bounds of several reactions should be set to zero according to supplementary material of 
constraints_2_zero = ['r66','r910','r1104','r239','r111','r106','r490','r791','r243','r252','r253','r307','r308','r404','r405','r1320','r639','r640','r641','r642','r649','r650','r651','r652','r645','r646','r643','r644','r653','r654','r655','r656','r534'] 
for rxn_id in constraints_2_zero:
    rxn = iLC915_model.reactions.get_by_id(rxn_id)
    rxn.bounds = (0,0)
# constraint formulate uptake because of not meaningful results
formulate_uptake = iLC915_model.reactions.get_by_id('r1143')
formulate_uptake.bounds = (0,0)
# important carbon sources (biotin, CO2, urea)
important_carbon = ['r1132','r1137','r1177']

# all carbon containing exchange reactions
carbon_exchange_rxns = ['r1122', 'r1123', 'r1124', 'r1126', 'r1127', 'r1129', 'r1130', 'r1131', 'r1132', 'r1134', 'r1135', 'r1137', 'r1138', 'r1139', 'r1140', 'r1141', 'r1144', 'r1145', 'r1146', 'r1147', 'r1148', 'r1149', 'r1151', 'r1152', 'r1153', 'r1154', 'r1155', 'r1156', 'r1157', 'r1158', 'r1161', 'r1162', 'r1163', 'r1165', 'r1167', 'r1168', 'r1170', 'r1172', 'r1173', 'r1174', 'r1175', 'r1176', 'r1177', 'r1178']

# set all carbon exchange reactions to 0 except important ones
for rxn_id in carbon_exchange_rxns:
    if rxn_id not in important_carbon:
        rxn = iLC915_model.reactions.get_by_id(rxn_id)
        rxn.bounds = (0,0)

# check definition of important rxns 
important_rxns = ['r1133', 'r1141', 'r1145', 'r1144', 'r1148', 'r1160', 'r1137', 'r1188']

if _verbose: 
    for rxn_id in important_rxns:
        rxn = iLC915_model.reactions.get_by_id(rxn_id)
        print(rxn.reactants, rxn.products, rxn.bounds)
        print(rxn_id, rxn.name, hf.formulaWithNames(rxn), 'rxn coefficient: ', rxn.get_coefficient(list(rxn.metabolites.keys())[0]))


# iterate over growth condition according to experimental data and store resulting growth rate
# rxns
model_name = 'iLC915'
biomass_reaction = 'r1133' # biomass
glu_ex_rxn = 'r1145' # glucose
fructose_ex_rxn = 'r1144' # fructose
glycerol_ex_rxn = 'r1148' # glycerol
o2_ex_rxn = 'r1160' # O2
co2_ex_rxn = 'r1137' # CO2
ethanol_ex_rxn = 'r1141' # ethanol
arabitol_ex_rxn = 'r1128' # Arabitol/arabinitol
maintenance = 'r1188' # maintenance
# set objective function
iLC915_model.objective = biomass_reaction

# growth results
iLC915_simulation_growth = []
# with CO2 constraint: 0.087; without CO2 constraint: 0.0898
iLC915_simulation_co2_prod = []

if (not _use_precomputed):
    with iLC915_model as model:
        for exp_idx in range(experiment_number):
            # set bounds for glucose and oxygen uptake
            glu_rxn = model.reactions.get_by_id(glu_ex_rxn) 
            glu_rxn.bounds = (glucose_uptake[exp_idx], glucose_uptake[exp_idx]) # rxn is defined as: "--> glucose"
            
            o2_rxn = model.reactions.get_by_id(o2_ex_rxn)
            o2_rxn.bounds = (oxygen_uptake[exp_idx], oxygen_uptake[exp_idx]) # rxn is defined as: "--> oxygen"

            # set maintenance
            maintenance_rxn = model.reactions.get_by_id(maintenance)
            # maintenance_rxn.bounds = (_ppa_maintenance,_ppa_maintenance)
            # maintenance_rxn.bounds = (1,1) # according to iLC915
            maintenance_rxn.bounds = (1.2434,1.2434) # according to iLC915 (1.2434 * 2.26 = 2.81)
            # maintenance_rxn.bounds = (0,0)

            # set ethanol production
            ethanol_rxn = model.reactions.get_by_id(ethanol_ex_rxn)
            if (ppa_ethanol_prod[exp_idx] == 1000): # reaction is defined as "ethanol -->"
                ethanol_rxn.bounds = (0,1000)
            else:
                ethanol_rxn.bounds = (ppa_ethanol_prod[exp_idx], ppa_ethanol_prod[exp_idx])
            

            # # set optional co2 uptake according to experimental data
            # co2_rxn = model.reactions.get_by_id(co2_ex_rxn)
            # co2_rxn.bounds = (ppa_co2_prod[exp_idx], ppa_co2_prod[exp_idx])

            # get solution
            solution = pfba(model)
            print(f'Growth rate on {glucose_uptake[exp_idx]} mmol/h glucose: {round(solution.fluxes[biomass_reaction],4)}')
            print(f'CO2 production on {glucose_uptake[exp_idx]} mmol/h glucose: {round(solution.fluxes[co2_ex_rxn],4)}')

            # store growth
            iLC915_simulation_growth.append(round(solution.fluxes[biomass_reaction],4))
            # store co2 production
            iLC915_simulation_co2_prod.append(round(solution.fluxes[co2_ex_rxn],4))

    # add results to value dicts
    _growth_conditions[model_name] = iLC915_simulation_growth
    _co2_productions[model_name] = iLC915_simulation_co2_prod
# _growth_conditions = {'experimental': [0.1, 0.1, 0.1], 'iLC915': [0.0898, 0.1089, 0.1302]}
# _co2_productions = {'experimental': [2.43, 2.55, 3.21], 'iLC915': [2.3224, 2.6005, 3.3094]}


# iMT1026: glucose results (ppa_glucose_consumption)
print('iMT1026v3 model: ')

# load model + set reactions
logging.getLogger("cobra").setLevel(logging.ERROR)
iMT1026v3_model = cobra.io.read_sbml_model(config['models']['ppa1026v3'])

# show interesting reactions
model_name = 'iMT1026v3'
biomass_reaction = 'growth' # default: 'growth' biomass
glu_ex_rxn = 'Ex_glc_D' # glucose
maintenance = 'ATPM' # maintenance (set to 2.9 in iMT1026v3)
glycerol_ex_rxn = 'Ex_glyc' # glycerol
o2_ex_rxn = 'Ex_o2' # O2
co2_ex_rxn = 'Ex_co2' # CO2
fructose_ex_rxn = 'Ex_fru' # fructose
ethanol_ex_rxn = 'Ex_etoh' # ethanol
arabitol_ex_rxn = 'Ex_abt_D' # Arabitol/arabinitol

interesting_rxns = ['growth', 'Ex_glc_D', 'ATPM', 'Ex_glyc', 'Ex_o2', 'Ex_co2', 'Ex_fru', 'Ex_etoh', 'Ex_abt_D']

# check definition of interesting rxns
if _verbose:
    for rxn_id in interesting_rxns:
        rxn = iMT1026v3_model.reactions.get_by_id(rxn_id)
        print(rxn.reactants, rxn.products, rxn.bounds)
        print(rxn_id, rxn.name, hf.formulaWithNames(rxn), 'rxn coefficient: ', rxn.get_coefficient(list(rxn.metabolites.keys())[0]), rxn.bounds)

# remove glyc reaction as carbon source
glyc_rxn = iMT1026v3_model.reactions.get_by_id('Ex_glyc')
glyc_rxn.bounds = (0, 1000)

try:
    iMT1026v3_model.summary() # not feasible => all good
    # raise error
    raise ValueError('Model behaves not as expected')
except:
    print('Model behaves as expected')

# iMT1026v3: iterate over growth condition according to experimental data and store resulting growth rate
# set objective function
iMT1026v3_model.objective = biomass_reaction

# growth results
iMT1026v3_simulation_growth = []
# with CO2 constraint: 0.087; without CO2 constraint: 0.0898
iMT1026v3_simulation_co2_prod = []
if (not _use_precomputed):
    with iMT1026v3_model as model:
        for exp_idx in range(experiment_number):
            # set bounds for glucose and oxygen uptake
            glu_rxn = model.reactions.get_by_id(glu_ex_rxn) 
            glu_rxn.bounds = (-glucose_uptake[exp_idx], -glucose_uptake[exp_idx]) # rxn is defined as: "--> glucose"
            
            o2_rxn = model.reactions.get_by_id(o2_ex_rxn)
            o2_rxn.bounds = (-oxygen_uptake[exp_idx], -oxygen_uptake[exp_idx]) # rxn is defined as: "--> oxygen"

            # set maintenance
            maintenance_rxn = model.reactions.get_by_id(maintenance)
            maintenance_rxn.bounds = (_ppa_maintenance,_ppa_maintenance)
            # maintenance_rxn.bounds = (0,0)
            # maintenance_rxn.bounds = (2.26,2.26)

            # set ethanol production
            ethanol_rxn = model.reactions.get_by_id(ethanol_ex_rxn)
            if (ppa_ethanol_prod[exp_idx] == 1000): # reaction is defined as "ethanol -->"
                ethanol_rxn.bounds = (0,1000)
            else:
                ethanol_rxn.bounds = (-ppa_ethanol_prod[exp_idx], -ppa_ethanol_prod[exp_idx])
            

            # # set optional co2 uptake according to experimental data
            # co2_rxn = model.reactions.get_by_id(co2_ex_rxn)
            # co2_rxn.bounds = (ppa_co2_prod[exp_idx], ppa_co2_prod[exp_idx])

            # get solution
            solution = pfba(model)
            print(f'Growth rate on {glucose_uptake[exp_idx]} mmol/h glucose: {round(solution.fluxes[biomass_reaction],4)}')
            print(f'CO2 production on {glucose_uptake[exp_idx]} mmol/h glucose: {round(solution.fluxes[co2_ex_rxn],4)}')

            # store growth
            iMT1026v3_simulation_growth.append(round(solution.fluxes[biomass_reaction],4))
            # store co2 production
            iMT1026v3_simulation_co2_prod.append(round(solution.fluxes[co2_ex_rxn],4))

    # add results to value dicts
    _growth_conditions[model_name] = iMT1026v3_simulation_growth
    _co2_productions[model_name] = iMT1026v3_simulation_co2_prod

if not _without_iMT1026v1:
    # iMT1026v1 (no meaningful results)
    print('iMT1026v1 model: ')

    logging.getLogger("cobra").setLevel(logging.ERROR)

    iMT1026v1_model = cobra.io.read_sbml_model(config['models']['ppa1026v1'])

    # reactions
    model_name = 'iMT1026v1'
    iMT1026v1_model.objective = 'growth' # default is BIOMASS
    biomass_reaction = 'growth'
    glu_ex_rxn = 'Ex_glc_D' # glucose
    maintenance = 'ATPM' # maintenance (set to 2.9 in iMT1026v3)
    glycerol_ex_rxn = 'Ex_glyc' # glycerol
    o2_ex_rxn = 'Ex_o2' # O2
    co2_ex_rxn = 'Ex_co2' # CO2
    fructose_ex_rxn = 'Ex_fru' # fructose
    ethanol_ex_rxn = 'Ex_etoh' # ethanol
    arabitol_ex_rxn = 'Ex_abt_D' # Arabitol/arabinitol

    # important_rxns = ['Ex_nh4', 'Ex_o2', 'Ex_pi', 'Ex_so4', 'Ex_co2', 'Ex_h2o', 'Ex_h', 'Ex_so3']
    important_rxns = ['Ex_nh4', 'Ex_o2', 'Ex_pi', 'Ex_so4', 'Ex_h2o', 'Ex_h', 'Ex_so3']

    for rxn in iMT1026v1_model.exchanges:
        if rxn.id not in important_rxns:
            rxn.bounds = (0,1000)

    # iMT1026v1: iterate over growth condition according to experimental data and store resulting growth rate
    if (not _use_precomputed):
        # set objective function
        iMT1026v1_model.objective = biomass_reaction

        # growth results
        iMT1026v1_simulation_growth = []
        # with CO2 constraint: 0.087; without CO2 constraint: 0.0898
        iMT1026v1_simulation_co2_prod = []

        with iMT1026v1_model as model:
            for exp_idx in range(experiment_number):
                # set bounds for glucose and oxygen uptake
                glu_rxn = model.reactions.get_by_id(glu_ex_rxn) 
                glu_rxn.bounds = (-glucose_uptake[exp_idx], -glucose_uptake[exp_idx]) # rxn is defined as: "--> glucose"
                
                o2_rxn = model.reactions.get_by_id(o2_ex_rxn)
                o2_rxn.bounds = (-oxygen_uptake[exp_idx], -oxygen_uptake[exp_idx]) # rxn is defined as: "--> oxygen"

                # set maintenance
                maintenance_rxn = model.reactions.get_by_id(maintenance)
                maintenance_rxn.bounds = (_ppa_maintenance,_ppa_maintenance)
                # maintenance_rxn.bounds = (0,0)

                # set ethanol production
                ethanol_rxn = model.reactions.get_by_id(ethanol_ex_rxn)
                if (ppa_ethanol_prod[exp_idx] == 1000): # reaction is defined as "ethanol -->"
                    ethanol_rxn.bounds = (0,1000)
                else:
                    ethanol_rxn.bounds = (-ppa_ethanol_prod[exp_idx], -ppa_ethanol_prod[exp_idx])
                

                # # set optional co2 uptake according to experimental data
                # co2_rxn = model.reactions.get_by_id(co2_ex_rxn)
                # co2_rxn.bounds = (ppa_co2_prod[exp_idx], ppa_co2_prod[exp_idx])

                # get solution
                solution = pfba(model)
                print(f'Growth rate on {glucose_uptake[exp_idx]} mmol/h glucose: {round(solution.fluxes[biomass_reaction],4)}')
                print(f'CO2 production on {glucose_uptake[exp_idx]} mmol/h glucose: {round(solution.fluxes[co2_ex_rxn],4)}')

                # store growth
                iMT1026v1_simulation_growth.append(round(solution.fluxes[biomass_reaction],4))
                # store co2 production
                iMT1026v1_simulation_co2_prod.append(round(solution.fluxes[co2_ex_rxn],4))

        # add results to value dicts
        _growth_conditions[model_name] = iMT1026v1_simulation_growth
        _co2_productions[model_name] = iMT1026v1_simulation_co2_prod


# # Results (precomputed):
# With iMT1026v1: 
if _use_precomputed: 
    _growth_conditions = {'experimental': [0.1, 0.1, 0.1],
    'iLC915': [0.0898, 0.1089, 0.1302],
    'iMT1026v3': [0.1028, 0.0894, 0.0738],
    'iMT1026v1': [-0.0, 0.0, 0.0016]}
    _co2_productions = {'experimental': [2.43, 2.55, 3.21],
    'iLC915': [2.3224, 2.6005, 3.3094],
    'iMT1026v3': [2.3999, 1.6526, 0.505],
    'iMT1026v1': [0.9452, 0.509, 0.1201]}

if _without_iMT1026v1:
    # remove iMT1026v1 from dicts
    try:
        del _growth_conditions['iMT1026v1']
        del _co2_productions['iMT1026v1']
    except:
        print('iMT1026v1 not in dicts')

# plot experimental vs simulated growth rates
title = f'Experimental vs simulated growth rates on glucose under different oxygen conditions for P. pastoris\n{_condition_name}'
y_label = 'Growth rate [1/h]'
outpath = ""
if _plotting:
    outpath = '../results/FBA_results/ppa_glucose_growth_comparison.png'
hf.plot_multiple_conditions(_growth_conditions, _condition_names, title, y_label, outpath, min_value = 0.05, max_value = 0.15, show_plot=True)


# plot experimental vs simulated co2 production
title = f'Experimental vs simulated CO2 production on glucose under different oxygen conditions for P. pastoris\n{_condition_name}'
y_label = 'CO2 production [mmol/gDW/h]'
outpath = ""
if _plotting:
    outpath = '../results/FBA_results/ppa_CO2_production_with_iMT1026v1.png'
    if _without_iMT1026v1:
        outpath = '../results/FBA_results/ppa_CO2_production_without_iMT1026v1.png'
    
hf.plot_multiple_conditions(_co2_productions, _condition_names, title, y_label, outpath, min_value = 0, max_value = 3.75, show_plot=True)