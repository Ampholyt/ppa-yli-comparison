# imports
import cobra
import numpy as np
import multiprocessing
import sys # append path
from cobra.flux_analysis import pfba # parsimonious FBA
import os 

sys.path.append('../scripts/')
import helperFunction as hf

config = hf.load_config()

# load model
model_name = 'yli21'
iYli21_model = cobra.io.read_sbml_model(config['models'][model_name])
glucose_rxn = 'R1070' # glucose
oxygen_rxn = 'R1287' # O2
biomass_rxn = 'biomass_C'

# shut glucose uptake off 
iYli21_model.reactions.get_by_id(glucose_rxn).bounds = (0, 0)

# growth rate on glucose: 0.24 h-1
# growth on glycerol: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3852309/ 0.32 h-1

_verbose = config['general']['verbose']
_snakemake = config['general']['snakemake']
min_glu = config['experiments']['PhPP']['min_glu']
# max_glu = 2.5 # max glucose uptake rate for pp according to literature (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320773/)
# max_glu = 15 # used for ecoli test model
# max_glu = 5 # just test for glucose 5 mmol (if set to this value 6.762 mmol o2 are required)
max_glu = config['experiments']['PhPP']['max_glu'] # just test for glucose 5 mmol (if set to this value 6.762 mmol o2 are required)
min_oxy = config['experiments']['PhPP']['min_oxy'] # default: 0.38
max_oxy = config['experiments']['PhPP']['max_oxy'] # default: 10
number_computations = config['experiments']['PhPP']['number_computations'] # default: 100


growth_rates = np.zeros((number_computations, number_computations))
glucose_range = np.linspace(min_glu, max_glu, number_computations)
oxygen_range = np.linspace(min_oxy, max_oxy, number_computations)


def comp_growth_rates(i):
    """Computes growth rates for a glucose value"""
    growth_rates_list = []
    for j in range(len(oxygen_range)):
        iYli21_model.reactions.get_by_id(glucose_rxn).bounds = (-glucose_range[i], -glucose_range[i])
        iYli21_model.reactions.get_by_id(oxygen_rxn).bounds = (-oxygen_range[j], oxygen_range[j])
        try:
            growth_rates_list.append(pfba(iYli21_model).fluxes[biomass_rxn])
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