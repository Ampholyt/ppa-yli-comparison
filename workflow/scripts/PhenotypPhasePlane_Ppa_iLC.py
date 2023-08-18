# imports
import cobra
import numpy as np
import multiprocessing
import sys # append path
from cobra.flux_analysis import pfba # parsimonious FBA
import os 


sys.path.append('../../scripts/')
import helperFunction as hf

config = hf.load_config()

# values: 
_verbose = config['general']['verbose']
_snakemake = config['general']['snakemake']
min_glu = config['experiments']['PhPP']['min_glu']
# max_glu = 2.5 # max glucose uptake rate for pp according to literature (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320773/)
# max_glu = 15 # used for ecoli test model
max_glu = config['experiments']['PhPP']['max_glu'] # just test for glucose 5 mmol (if set to this value 9.5 mmol o2 are required)
min_oxy = config['experiments']['PhPP']['min_oxy'] # default: 0.38
max_oxy = config['experiments']['PhPP']['max_oxy'] # default: 10
number_computations = config['experiments']['PhPP']['number_computations'] # default: 100


# load model
model_name = 'ppaiLC915'
iLC915_model = cobra.io.read_sbml_model(config['models'][model_name])
oxygen_rxn = 'r1160' # O2
glucose_rxn = 'r1145' # glucose
maintenance = 'r1188' # maintenance
biomass_rxn = 'r1339' # growth


# prepare model
print(f'prepare {model_name} model')
# prepare iLC915_model: needed because of internal cycles and active reactions
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


maintenance_rxn = iLC915_model.reactions.get_by_id(maintenance)

maintenance_rxn.bounds = (1.2434,1.2434) # according to iLC915 (1.2434 * 2.26 = 2.81)

# set objective: 
iLC915_model.objective = biomass_rxn


growth_rates = np.zeros((number_computations, number_computations))
glucose_range = np.linspace(min_glu, max_glu, number_computations)
oxygen_range = np.linspace(min_oxy, max_oxy, number_computations)


def comp_growth_rates(i):
    """Computes growth rates for a glucose value"""
    growth_rates_list = []
    for j in range(len(oxygen_range)):
        iLC915_model.reactions.get_by_id(glucose_rxn).bounds = (glucose_range[i], glucose_range[i]) # --> alpha-D-Glucose_C6H12O6
        iLC915_model.reactions.get_by_id(oxygen_rxn).bounds = (-oxygen_range[j], oxygen_range[j]) # # --> Oxygen_O2
        try:
            growth_rates_list.append(pfba(iLC915_model).fluxes[biomass_rxn]) 
        except:
            growth_rates_list.append(0)
        if _verbose:
            print('glucose: {}, oxygen: {}, growth rate: {}'.format(glucose_range[i], oxygen_range[j], growth_rates_list[j]))
    return growth_rates_list


def main():
    print('computing ...')
    # create pool
    pool = multiprocessing.Pool(processes = 8)
    returns = pool.map(comp_growth_rates, range(len(glucose_range)))

    print(f'computing for {model_name} model finished')

    # store data
    condition_name = f'glu_{max_glu}_oxy_{max_oxy}'
    outdir = 'results/PhPP/'
    file_name = "{}x{}_{}_{}_growth_rates".format(number_computations, number_computations, condition_name, model_name)
    outfile = outdir + file_name
    if _snakemake:
        outfile = snakemake.output[0]
    np.save(outfile, returns)
    if _verbose:
        print(returns)


if __name__ == '__main__':
    main()