# imports
import sys
import os
import pandas as pd
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis.parsimonious import pfba

import helperFunction as hf # helper functions

# This file:
# a) calculates parsimonious FBA and adds growth rates into NumbersDF
# b) add or creates a flux dataframe with a certain name (comparison of different flux conditions) but it is extendable to prepare a flux frame for different models (need to connect the reactions)
# 
# 2 Models:
# a) Given: model, outputDir: It takes the default values (glu, 'R1070' and (-2.43, 1000))
# expecting: 
# python FluxAnalysis.py <path to model> <output directory> <conditionName> <exchange reactions> <lower bounds> <upper bounds>
# Exp: 
# Version1) python3 FluxAnalysis.py  /Users/ampholyt/Coding/BSEP22/Data/Models/iYli21_v1.xml /Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/
# Version2) python3 FluxAnalysis.py  /Users/ampholyt/Coding/BSEP22/Data/Models/iYli21_v1.xml /Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/ gly R1141 -2.43 1000
####### default cases: 
_condition = 'glu' # glucose
_exchangeIds = ['R1070']
_exchangeBounds = [(-2.43,1000)]
_sep=';'
######

# Interact with command line arguments
n = len(sys.argv)
_verbose = True
if (n != 3):
    if (n != 7):
        raise Exception('Wrong number of arguments detected.')
    _condition = sys.argv[3]
    _exchangeIds = sys.argv[4].split(',')
    _lowerBounds = sys.argv[5].split(',')
    _upperBounds = sys.argv[6].split(',')
    # Allow modifying multiple exchange reactions
    _exchangeBounds = [(float(_lowerBounds[idx]),float(_upperBounds[idx])) for idx in range(len(_lowerBounds))]
    print(_exchangeBounds)

# paths
_modelPath = sys.argv[1]
_outputDir = sys.argv[2]

_modelPrefix = _modelPath.split('/')[-1].split('.')[0]
_fluxDf = '/Fluxes/{modelPrefix}_fluxDf.csv'.format(modelPrefix=_modelPrefix)
_numbersDf = '/Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/ModelOverview/{modelPrefix}_modelNumbers.csv'.format(modelPrefix=_modelPrefix)
_defaultFile = '/Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/ModelOverview/'+ _modelPrefix.split('.')[0] + '_exchangeReactionBounds.csv'
if _outputDir[-1] != '/':
    _outputDir += '/'

print(f"Assuming that the fluxDf should be stored in {_outputDir + '/Fluxes/'} and ModelOverview.py ran already.")
print('Flux dataframe: ', _fluxDf)
print('Numbers dataframe: ', _numbersDf)


def pFBAofCondition(model, defaultFile, exchangeIds, bounds):
    '''Generartes model and pFBA result based on exchangeIds and bounds'''
    # default exchange condition
    hf.exchangeToDefault(model, defaultFile) # without exchange of carbon
    if (len(exchangeIds) != len(bounds)):
        raise Exception('ExchangeIds and Bounds need to have the same length.')
    # Set bounds of exchange reactions (to -2.43 mmol/h based on Guo et al. 2022)
    for idx in range(len(exchangeIds)):
        exchangeReac = model.reactions.get_by_id(exchangeIds[idx])
        exchangeReac.bounds = bounds[idx]
    return pfba(model)

def addGrowthRate(model, defaultFile, conditionName, numberDf, exchangeIdList, boundList):
    # generate growth rate
    pfba_solution = pFBAofCondition(model, defaultFile, exchangeIdList, boundList)
    growthRate = pfba_solution.fluxes.biomass_C
    # load numbers dataframe and add growth rate on condition
    numbersDf = pd.read_csv(numberDf, sep=_sep)
    numbersDf[f'{conditionName}Growth'] = [growthRate]
    numbersDf.to_csv(numberDf, sep=_sep)

def createFluxFrame(solution, pathwayDf, outDir, condition):
    '''Creates a new flux Frame based on the available data'''
    fluxSer = solution.fluxes
    fluxDf = pd.DataFrame({'reactionId': fluxSer.index, f'{condition}Flux': fluxSer})
    # merge flux frame with pathway frame
    fluxPathwayDf = fluxDf.merge(pathwayDf, on='reactionId', how='left')
    fluxPathwayDf = fluxPathwayDf[['reactionId', 'pathway', f'{condition}Flux']]
    fluxPathwayDf.to_csv(outDir, sep=_sep, index = False)
    return fluxPathwayDf

def addFluxToFrame(solution, outDir, condition):
    '''Adds calculated fluxes into the existing flux frame'''
    fluxDf = pd.read_csv(outDir, sep=_sep)
    fluxSer = solution.fluxes
    # needed for the additional series information
    fluxDf.index = fluxDf['reactionId']
    fluxDf[f'{condition}Flux'] = fluxSer
    fluxDf.to_csv(outDir, sep=_sep, index = False)

def createPathwayFrame(model):
    '''Prepare dataframe from pathway names and reaction ids'''
    pathways = model.groups
    reactionIds = []
    pathwayNames = []
    for pathway in pathways:
        for reaction in pathway.members:
            pathwayNames.append(pathway.name)
            reactionIds.append(reaction.id)
    pathwayDf = pd.DataFrame({'pathway':pathwayNames, 'reactionId':reactionIds})
    return pathwayDf

def extendFluxFrame(model, defaultFile, fluxDfDir, conditionName, exchangeIds, bounds):
    '''Creates or Adds flux data frame based on condition (only one exchangeReaction is possible)'''
    pfba_solution = pFBAofCondition(model, defaultFile, exchangeIds, bounds)
    if os.path.exists(fluxDfDir):
        addFluxToFrame(pfba_solution, fluxDfDir, conditionName)
    else:
        # create pathway frame
        pathwayDf = createPathwayFrame(model)
        createFluxFrame(pfba_solution, pathwayDf, fluxDfDir, conditionName)

def main():
    # read model
    model = read_sbml_model(_modelPath)
    # growth rates:
    addGrowthRate(model, _defaultFile, _condition, _numbersDf, _exchangeIds, _exchangeBounds)
    # create or add fluxes to fluxes dataframe
    extendFluxFrame(model, _defaultFile, _outputDir + _fluxDf, _condition, _exchangeIds, _exchangeBounds)
    print(sys.argv[0], 'finished.')

if __name__ == '__main__':
    main()