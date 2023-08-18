import cobra
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
import pandas as pd

model = cobra.io.read_sbml_model('iYli21_v1.sbml')
_modelPrefix = 'iYli21_v1'
_condition = 'glu'
_defaultFile = _modelPrefix.split('.')[0] + '_exchangeReactionBounds.csv'


def setReactionBounds(model, boundFile):
    '''Gets a file (expecting at least 2 columns: "reactionId" and "bounds") 
    with bound information for the exchange reactions and sets the corresponding reaction bounds.'''
    for i, reactionId in boundFile['reactionId'].items():
        # get default reaction bounds (eval (tuple))
        reacBounds = eval(boundFile['bounds'][i])
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
    '''Reading default file and setting model to default'''
    print(f'Set reaction bounds from model {model} to default')
    # Read file line by line
    try:
        defaultDf = pd.read_csv(defaultFile)
    except:
        raise Exception('Default values are not defined:\n' + defaultFile + '\nNo such file or directory.')
    # set every exchange reaction to zero: assuming echange in reaction name
    setExchangeToZero(model)
    # set reaction bounds
    setReactionBounds(model,defaultDf)

# set condition to glucose
exchangeToDefault(model, _defaultFile) # without exchange of carbon

# Set glucose reaction to -2.43 mmol/h (based on Guo et al. 2022)
dGlu = model.reactions.get_by_id('R1070')
dGlu.bounds = (-2.43, 1000)


# 1. single gene deletion
deletion_results = single_gene_deletion(model)
deletion_results.to_csv(f'{_condition}SDeletion_{_modelPrefix}.csv', header=True)

# 2. single reaction deletion
rea_deletion_results = single_reaction_deletion(model)
rea_deletion_results.to_csv(f'{_condition}RDeletion_{_modelPrefix}.csv', header=True)

# 3. double gene deletion
doubleG_deletion_results = double_gene_deletion(model)
doubleG_deletion_results.to_csv(f'{_condition}DouSDeletion_{_modelPrefix}.csv', header=True)

# 4. double reaction deletion
doubleR_deletion_results = double_reaction_deletion(model)
doubleR_deletion_results.to_csv(f'{_condition}DouRDeletion_{_modelPrefix}.csv', sep='\t', header=True)