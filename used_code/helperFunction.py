# helper functions:
import cobra
from configparser import ConfigParser, ExtendedInterpolation # config file parsing
import os
import pandas as pd
import numpy as np # plot multiple conditions
import matplotlib.pyplot as plt # plot multiple conditions
import yaml # config file parsing

_sep = ';'


# print functions or extensions of implemented features
def reactionInformationPrint(reactions, id_name_same = False):
    '''Prints information about the reaction in a readable format'''

    for reaction in reactions:

        if id_name_same: # print name of first metabolite
            print(f'{list(reaction.metabolites)[0].name}, (id: {reaction.id}) is active. (Bounds: {reaction.bounds})')
            continue

        print('{name}, (id: {id}) is active. (Bounds: {bounds}); Reaction: {reaction}'.format(name = reaction.name, id = reaction.id, bounds=reaction.bounds, reaction = reaction.reaction))

def mediumAnalysis(model, id_name_same = False):
    '''Prints the current active media reactions in a readable format'''
    newModel = model.copy()
    mediaReactions = []
    print('\n----- media analysis -----')
    for key,value in newModel.medium.items():
        # search for keys and return names
        mediaReactions.append(newModel.reactions.get_by_id(key))
    reactionInformationPrint(mediaReactions, id_name_same)
    print('--------------------------')

def formulaWithNames(reaction):
    '''Returns the reaction formula with metabolite names instead of ids'''
    formula = reaction.reaction
    for metabolite in reaction.metabolites:
        formula = formula.replace(metabolite.id, metabolite.name)
    return formula

## print the metabolites dict nicely
def metabolites_print(model, rxn_id):
    """Displays the metabolites and their coefficients human readable in the console."""
    rxn = model.reactions.get_by_id(rxn_id)
    print(f'Metabolites for the reaction: {rxn.name} (id: {rxn.id})')
    for metabolite, coefficient in rxn.metabolites.items():
        print(f'{metabolite.name} ({metabolite.id}): {coefficient}')
    return True

# New Features
## Set Exchange to default 
def setReactionBounds(model, boundFile):
    '''Gets a file (expecting at least 2 columns: "reactionId" and "bounds") 
    with bound information for the exchange reactions and sets the corresponding reaction bounds.'''
    for i, reactionId in boundFile['reactionId'].items():
        # get default reaction bounds (eval (tuple))
        reacBounds = eval(boundFile['reactionBounds'][i])
        # get model reaction with reactionid
        modelReaction = model.reactions.get_by_id(reactionId)
        modelReaction.bounds = reacBounds

def setExchangeToZero(model):
    for reaction in model.reactions:
        # find reactions with exchange in reaction name
        if 'exchange' in reaction.name:
            # set bounds to zero
            reaction.bounds = (0.0,0.0)

def exchangeToDefault(model, defaultFile): 
    # TODO: change depending on file make object of bounds also possible
    '''Reading default file and setting model to default'''
    print(f'Set reaction bounds from model {model} to default')
    # Read file line by line
    try:
        defaultDf = pd.read_csv(defaultFile, sep=_sep)
    except:
        raise Exception('Default values are not defined:\n' + defaultFile + '\nNo such file or directory.')
    # set every exchange reaction to zero: assuming echange in reaction name
    setExchangeToZero(model)
    # set reaction bounds
    setReactionBounds(model,defaultDf)

## Set conditions to Model
### Used in FluxAnalysis and essential Gene analysis
def conditionToModel(model, defaultFile, exchangeIds, bounds):
    '''Sets given condition to model (exchangeIds and bounds)'''
    # default exchange condition
    defModel = model.copy()
    exchangeToDefault(defModel, defaultFile) # without exchange of carbon
    if (len(exchangeIds) != len(bounds)):
        raise Exception('ExchangeIds and Bounds need to have the same length.')
    # Set bounds of exchange reactions (default: to -2.43 mmol/h based on Guo et al. 2022)
    for idx in range(len(exchangeIds)):
        exchangeReac = defModel.reactions.get_by_id(exchangeIds[idx])
        exchangeReac.bounds = bounds[idx]
    return defModel

def setConditionToModel(model, condDf, condition, defaultFile):
    '''Sets the given condition to the model. 
        Parameters: model, condition dataframe, condition name, default file'''
    condDf = condDf[condDf['conditionName'] == condition]
    reactionId = eval(condDf.exchangeReactions.values[0])
    bounds = eval(condDf.bounds.values[0]) # converts string to tuple
    return conditionToModel(model, defaultFile, reactionId, bounds)

def read_config(configPath):
    config = ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configPath)
    return config

def getConditionDf(configPath):
    '''Load condition df based on condition.csv file given in config file'''
    # 2.1. find the uptake reaction for polyol (glycerol), sugar (glucose), fatty acid (oleic acid), alkane (hexadecane) and glycerolipids (triolein) and (tributyrin)
    config = read_config(configPath)
    sep = config['DEFAULT']['separator']
    # read csv file
    conditionDf = pd.read_csv(config['DEFAULT']['condition_path'], sep=sep)
    return conditionDf

# Interacting with metabolite, gene or reaction lists
def getIds(elementList):
    '''Returns the ids of an element list (list of genes, reactions or metabolites'''
    return [elem.id for elem in elementList]

def getNames(elementList):
    '''Returns the names of the elemts in the given list (list of genes, reactions or metabolites'''
    return [elem.name for elem in elementList]

# Find all reactions which are associated with a given gene (or list of genes)
def findReactionsForGenes(model, gene_list, blocked_rxn_list = []):
    '''Returns a list of reactions which are associated with the given genes and not in blocked_rxn_list'''
    # find all reactions which are associated with the given genes
    essential_reactions = set()
    for gene_id in gene_list:
        for reaction in model.genes.get_by_id(gene_id).reactions:
            # check if the reaction is in the blocked reaction list
            if reaction.id not in blocked_rxn_list: 
                essential_reactions.add(reaction.id)
    return list(essential_reactions)

def model_specific_manipulations(model_name, ex_rxn):
    """Takes model name and exchange reaction name and returns the exchange metabolite name."""
    if model_name == 'ppa1026v3':
        exchange_name = ex_rxn.name.lower().split(' exchange')[0]
        # normalize names: remove ', yeast-specific'; "3',5'-" => 3,5; '1,3' => 1-3
        if ', yeast-specific' in exchange_name:
            exchange_name = exchange_name.split(', yeast-specific')[0]
        exchange_name = exchange_name.replace(' ', '_').replace('-', '_').replace("'", '').replace(',', '-')
    
    elif model_name == 'yli21':
        exchange_name = ex_rxn.name
        # 4 cases where ' exchange' not given EXC_OUT_m1803, EXC_OUT_m1826, EXC_OUT_m1824 (isocitrate_C6H8O7, erythritol_, D-mannitol_) + " transport" => remove "_x"
        if 'EXC_OUT_' in exchange_name:
            exchange_name = list(ex_rxn.metabolites)[0].name.split('_')[0]
        elif ' transport' in exchange_name:
            # remove ' transport'
            exchange_name = exchange_name.split(' transport')[0]
        else:
            # remove ' exchange'
            exchange_name = ex_rxn.name.split(' exchange')[0]
        
        # lower; ',' => -; "'" => ''; ' ' => '_'
        exchange_name = exchange_name.lower().replace(' ', '_').replace('-', '_').replace("'", '').replace(',', '-')
        
    elif model_name == 'yli4':
        # one case with transport
        exchange_name = ex_rxn.name
        # 4 cases where ' exchange' not given EXC_OUT_m1803, EXC_OUT_m1826, EXC_OUT_m1824 (isocitrate_C6H8O7, erythritol_, D-mannitol_) + " transport" => remove "_x"
        if 'EXC_OUT_' in exchange_name:
            exchange_name = list(ex_rxn.metabolites)[0].name.split('_')[0]
        elif ' transport' in exchange_name:
            # remove ' transport'
            exchange_name = exchange_name.split(' transport')[0]
        else:
            # remove ' exchange'
            exchange_name = exchange_name.split(' exchange')[0]
        
        # lower; ',' => -; "'" => ''; ' ' => '_'
        exchange_name = exchange_name.lower().replace(' ', '_').replace('-', '_').replace("'", '').replace(',', '-')
    return exchange_name

def exchange_metabolite_table(exchangeRxn, model_name, output_path, sep, verbose=False):
    """Prepares the exchange metabolite table for a given list of exchange reactions. Returns a pandas dataframe with the exchange metabolite name and id."""
    metabolite_names = []
    metabolite_ids = []
    unnormalized_names = []
    for ex_rxn in models[model].exchanges:
        try:
            list(ex_rxn.metabolites)[1]
            raise Exception('more than one metabolite found, check that models exchange reactions in more detail')
        except:
            exchange_name = model_specific_manipulations(model_name, ex_rxn)
            metabolite_names.append(exchange_name)
            unnormalized_names.append(ex_rxn.name)
            metabolite_ids.append(list(ex_rxn.metabolites)[0].id)      
        
    if (len(metabolite_names) != len(set(metabolite_ids))):
        raise Exception('metabolite names and ids are not unique')
    # convert to two lists to dataframe with two columns (exchange_metabolite_name, exchange_metabolite_id)
    exchange_metabolites_df = pd.DataFrame({'exchange_metabolite_name': metabolite_names, 'exchange_metabolite_id': metabolite_ids, 'unnormalized_name':  unnormalized_names})
    # save to csv
    out_file_path = os.path.join(output_path,f'{model}_exchange_metabolites.csv')
    # make sure the output path exists
    os.makedirs(os.path.dirname(out_file_path), exist_ok=True)
    exchange_metabolites_df.to_csv(out_file_path, index=False, sep=sep)
    if verbose:
        print(f'Exchange metabolite table of {model_name} was written to: {out_file_path}')
    return exchange_metabolites_df

# set model conditions
def prepare_iLC915_model(model):
    """Gets the model and prepares it for the simulation.
    - Remove every unnecessary reaction
    - Setup the minimal medium
    - Setup Glucose uptake
    """
    # remove all reactions
    for rxn in model.exchanges:
        rxn.bounds = (0,0)
    
    for rxn in model.boundary:
        rxn.bounds = (0,0)

    # setup the important reactions (H2O, Sulfate, Fe2+, NH3, Orthophosphate)
    min_medium_boundary = ['r1150', 'r1169', 'r1142', 'r1159', 'r1164']
    
    # limit as given in the experimental data
    oxygen = ['r1160']

    # replace (0,0) bounds with (0,1000) bounds to allow growth
    biomass = ['r1133']

    # glucose as carbon source 
    glucose = ['r1145'] 
    glu_uptake = 1 # experimental: 1 and 1.72 mmol/h
    # allow growth: 
    for rxn in model.exchanges:
        if rxn.id in min_medium_boundary:
            rxn.bounds = (-1000, 1000)
        elif rxn.id in oxygen:
            rxn.bounds = (-1000, 2.35)
        elif rxn.id in biomass:
            rxn.bounds = (0,1000)
        elif rxn.id in glucose:
            rxn.bounds = (0,glu_uptake)

def getReactionFluxes(model, reactions, biomass_rxn, growth_condition, maintenance_rxn, maintenance, verbose=False):
    """Computes FBA solution of a given growth condition and returns the fluxes of the reactions of interest"""
    current_fluxes = {}

    # compute FBA solution
    # set objective function
    model.objective = biomass_rxn
    # set growth condition
    for rxn_id, flux in growth_condition.items():
        model.reactions.get_by_id(rxn_id).lower_bound = -flux
    # set maintenance
    model.reactions.get_by_id(maintenance_rxn).bounds = maintenance
    # pfba solution
    pfba_solution = cobra.flux_analysis.pfba(model)
    if verbose:
        print(f'The objective value of the {model} is: {pfba_solution.fluxes[biomass_rxn]}')
    # collect pathway fluxes
    for ec_number, rxn_id in reactions.items():
        current_fluxes[ec_number] = (pfba_solution.fluxes[rxn_id],rxn_id)
    return current_fluxes

def plot_multiple_conditions(value_dict, conditions, title, y_label, outpath, min_value = 0, max_value = 0.3, show_plot = True):
    '''Plot the for each condition of multiple models the given values (in a dict of value lists) + list of conditions
    e.g. value_dict = {'iYali4_model': [0.0036, 0.0045],...}
    @params: conditions: list of conditions
    @trick: (is extendable) is able to plot conditions for different models'''
    # get min and max values of dict of lists
    # min_value = max(0, min([sublist[-1] for sublist in list(value_dict.values())]))
    # max_value = max([sublist[-1] for sublist in list(value_dict.values())])

    x = np.arange(len(conditions))  # the label locations
    width = 0.1  # the width of the bars
    multiplier = 0
    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in value_dict.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute.split('_')[0])
        # ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.set_xticks(x + width, conditions)
    ax.legend(loc='upper left', ncol=3)
    ax.set_ylim(min_value, max_value)
    if outpath != "": 
        plt.savefig(outpath)
    if show_plot:
        plt.show()

def numberDf(model, modelPrefix, outputDir, sep, createOutput, verbose):
    '''Check numbers of genes, metabolites and reactions and safes them to dataframe'''
    # investigate number and print numbers
    modelName = model.name
    reacNum = len(model.reactions)
    metaboNum = len(model.metabolites)
    geneNum = len(model.genes)
    if verbose:
        print('The {name} has {reacNum} reactions, {metaboNum} metabolites and {geneNum} genes'.format(name=modelName, reacNum=reacNum, metaboNum=metaboNum, geneNum=geneNum))
    # generate dataframe
    tupleList = [(modelPrefix,geneNum,metaboNum,reacNum)]
    modelNumbers = pd.DataFrame(tupleList, columns=['modelName','#Genes','#Metabolites','#Reactions'])
    if createOutput:
        if verbose:
            print('The model numbers will be saved to a .csv file.', "Path: ", outputDir+modelPrefix+'_modelNumbers.csv')
        modelNumbers.to_csv(outputDir + modelPrefix + '_modelNumbers.csv', index=False, sep=sep)
    return modelNumbers

def load_config(config_name='model_config'):
    config_path = f'../config/{config_name}.yaml'

    # load config
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def get_all_models():
    """Returns a list of all models in use (also S. cerevisiae and E. coli in order to compare them)"""
    config = load_config()
    models = {
        'yli647_corr': config['models']['yli647_corr'], 
        'yli647_uncorr': config['models']['yli647_uncorr'], 
        'iYali4_corr': config['models']['yli4_corr'], 
        'iYali4_uncorr': config['models']['yli4_uncorr'], 
        'iYli_2.0_corr': config['models']['yli2.0_corr'], 
        #'iYli_2.0_uncorr': config['models']['yli2.0_uncorr'], 
        'iMK735_corr': config['models']['yliMK735_corr'], 
        'iMK735_uncorr': config['models']['yliMK735_uncorr'], 
        'iNL895_corr': config['models']['yliNL895_corr'],
        'iNL895_uncorr': config['models']['yliNL895_uncorr'],
        'iYli21': config['models']['yli21'],
        'iYL619_PCP': config['models']['iYL619_PCP'],
        'PpaMBEL1254': config['models']['ppaMBEL1254'],
        'iMT1026v3': config['models']['ppa1026v3'],
        'iMT1026Chan2017': config['models']['ppa1026Chan'],
        'iLC915': config['models']['ppaiLC915'],
        'iAF1260': config['models']['ecoli'],
        'iMM904': config['models']['scere'],
    }
    return models

def get_loadable_models():
    config = load_config()

    loadable_models = {
        'yli647_corr': config['models']['yli647_corr'],
        'iYali4_corr': config['models']['yli4_corr'], 
        'iYali4_uncorr': config['models']['yli4_uncorr'], 
        'iYli_2.0_corr': config['models']['yli2.0_corr'], 
        'iMK735_corr': config['models']['yliMK735_corr'], 
        'iMK735_uncorr': config['models']['yliMK735_uncorr'], 
        'iNL895_corr': config['models']['yliNL895_corr'],
        'iNL895_uncorr': config['models']['yliNL895_uncorr'],
        'iYli21': config['models']['yli21'],
        'PpaMBEL1254': config['models']['ppaMBEL1254'],
        'iMT1026v3': config['models']['ppa1026v3'],
        'iLC915': config['models']['ppaiLC915'],
    }
    return loadable_models

def get_models():
    """Returns a list of the models from Yli and Ppa in use"""
    config = load_config()
    models = {
        'yli647_corr': config['models']['yli647_corr'], 
        'iYali4_corr': config['models']['yli4_corr'], 
        'iYli_2.0_corr': config['models']['yli2.0_corr'], 
        'iMK735_corr': config['models']['yliMK735_corr'], 
        'iNL895_corr': config['models']['yliNL895_corr'],
        'iYli21': config['models']['yli21'],
        'iMT1026v3': config['models']['ppa1026v3'],
        'iMT1026Chan2017': config['models']['ppa1026Chan'],
        'iLC915': config['models']['ppaiLC915'],
    }
    return models