## This skrip loads the exchange reactions from all the models and compares them 
#### Which are exchange reactions in all models, What is the id of glucose and fructose for each model
#### As a sidenote: return list of reactions

# imports
import cobra
import yaml
import os
import pandas

import helperFunction as hf


config_name = 'model_config'
# config
config_path = f'../config/{config_name}.yaml'

# outpath
outpath = '../results/exchange_metabolites' # without '/' at the end

# load config
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# load models
model_names = list(config['models'].keys())
_sep = config['seperator']
yeast_models = ['ppa1026v3', 'yli21', 'yli4']
models = {}

for model, path in config['models'].items():
    models[model] = cobra.io.read_sbml_model(path)

def check_number_of_exchange_reactions(models_dict, yeast_models):
    """Print number of exchange reactions for each model (iMT1026v3, iYli21, iYli4)"""
    iMT_exchange = models[yeast_models[0]].exchanges
    iYli21_exchange = models[yeast_models[1]].exchanges
    iYli4_exchange = models[yeast_models[2]].exchanges
    print('iMT1026v3: ', len(iMT_exchange))
    print('iYli21: ', len(iYli21_exchange))
    print('iYali4: ', len(iYli4_exchange))

# main function
def main():
    # check number of exchange reactions
    check_number_of_exchange_reactions(models, yeast_models)

    # iterate over models and transform exchange reactions into a table
    for model in yeast_models:
        exchange_table = hf.exchange_metabolite_table(models[model].exchanges, model, outpath, _sep)
    
    # check how many exchange reactions do the models have in common
    ## load exchange reactions dataframes
    iYli21_exchange_df = pandas.read_csv(f'{outpath}/yli21_exchange_metabolites.csv', sep=_sep)
    iMT1026_exchange_df = pandas.read_csv(f'{outpath}/ppa1026v3_exchange_metabolites.csv', sep=_sep)
    iYali4_exchange_df = pandas.read_csv(f'{outpath}/yli4_exchange_metabolites.csv', sep=_sep)

    # get sets of exchange reactions
    iYli21_exchange_set = set(iYli21_exchange_df['exchange_metabolite_name'])
    iMT1026_exchange_set = set(iMT1026_exchange_df['exchange_metabolite_name'])
    iYali4_exchange_set = set(iYali4_exchange_df['exchange_metabolite_name'])

    # intersection of all three sets
    iYli21_iMT1026_iYali4_common = iYli21_exchange_set.intersection(iMT1026_exchange_set, iYali4_exchange_set)
    print('iYli21_iMT1026_iYali4_common: ', len(iYli21_iMT1026_iYali4_common)) # 95

    # find the ids of the metabolites in common
    iYli21_common_ids = []
    iYali4_common_ids = []
    iMT1026_common_ids = []
    unnormalized_names = []
    for metabolite_name in iYli21_iMT1026_iYali4_common:
        iYli21_common_ids.append(iYli21_exchange_df[iYli21_exchange_df['exchange_metabolite_name'] == metabolite_name]['exchange_metabolite_id'].values[0])
        iYali4_common_ids.append(iYali4_exchange_df[iYali4_exchange_df['exchange_metabolite_name'] == metabolite_name]['exchange_metabolite_id'].values[0])
        iMT1026_common_ids.append(iMT1026_exchange_df[iMT1026_exchange_df['exchange_metabolite_name'] == metabolite_name]['exchange_metabolite_id'].values[0])
        unnormalized_names.append(iMT1026_exchange_df[iMT1026_exchange_df['exchange_metabolite_name'] == metabolite_name]['unnormalized_name'].values[0])

    # create dataframe with metabolite names and ids
    iYli21_iMT1026_iYali4_common_df = pandas.DataFrame({'exchange_metabolite_name': list(iYli21_iMT1026_iYali4_common), 'exchange_metabolite_id_iYli21': iYli21_common_ids, 'exchange_metabolite_id_iYali4': iYali4_common_ids, 'exchange_metabolite_id_iMT1026': iMT1026_common_ids, 'exchange_reaction_name': unnormalized_names})
    # save to csv
    out_file_path = os.path.join(outpath,f'common_exchange_metabolites_iYli21_ids.csv')
    # make sure the output path exists
    os.makedirs(os.path.dirname(out_file_path), exist_ok=True)
    # store the dataframe of common exchange reactions
    iYli21_iMT1026_iYali4_common_df.to_csv(out_file_path, index=False, sep=_sep)

if __name__ == '__main__':
    main()