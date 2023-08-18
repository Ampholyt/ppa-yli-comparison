# This skript will load the flux frame file (flux distributions for all reactions) for all conditions of one Model.
# It will visualize the results of filtering based on user specified values. (Heatmaps)
# expecting fluxVisluaization.py model=<model>(req) prefix=<modelPrefix> cond=<condition>(req) flux=<fluxFrame>(req) filter=<filterThreshold> <outputDir>(req)

# imports:
import cobra
from cobra.io import read_sbml_model
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# default: 
_modelPath = '/Users/ampholyt/Coding/BSEP22/Data/Models/iYli21_v1.xml'
_outputDir = '/Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/'
_heatDir = '/Heatmaps/'
_modelPrefix = 'iYli21_v1'
_condition = 'glu'
_fluxDf = _outputDir + f'/Fluxes/{_condition}_fluxDf.csv'

print('Warning: This script assumes so called interesting pathways which are model specific')

# setup:
# sys.argv = ['fluxVisualization.py', 'model=iYli21_v1', 'prefix=iYli21_v1', 'cond=glu', 'flux=/Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/Fluxes/glu_fluxDf.csv', 'filter=0.1', 'outdir=/Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/']

# argCounter = 0
# for arg in sys.argv:
#     if 'model=' in arg:
#         _modelPath = arg.split('=')[1]
#         argCounter += 1
#     elif 'prefix=' in arg:
#         _modelPrefix = arg.split('=')[1]
#     elif 'cond=' in arg:
#         _condition = arg.split('=')[1]
#         argCounter += 1
#     elif 'flux=' in arg:
#         _fluxDf = arg.split('=')[1]
#         argCounter += 1
#     elif 'filter=' in arg:
#         _filterThreshold = float(arg.split('=')[1])
#     elif 'outdir=' in arg:
#         _outputDir = arg.split('=')[1]
#         _heatDir = _outputDir + '/Heatmaps/'
#         argCounter += 1
    
# if argCounter < 4:
#     raise Exception('Not enough arguments provided. Please provide at least model, condition, fluxDf and outputDir.')

# _fluxDf = '/Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/Fluxes/glu_fluxDf.csv'


# load helper functions

def filteredFlux(fluxDf, filterThreshold, outputDir, modelPrefix, condition):
    # heatmap of filtered flux
    filteredFlux = fluxDf.copy()

    # calculate max, min and 95% percentil of fluxes
    quantileVal = fluxDf[f'{condition}Flux'].quantile(0.95)
    maxFlux = fluxDf[f'{condition}Flux'].max()
    minFlux = fluxDf[f'{condition}Flux'].min()
    print(f'95% quantile of fluxes: {quantileVal}\nMax flux: {maxFlux}\nMin flux: {minFlux}')

    # filter flux values for absolute fluxes > _filterThreshold
    filteredFlux['absFlux'] = filteredFlux[f'{condition}Flux'].abs()
    filteredFlux = filteredFlux[filteredFlux['absFlux'] > filterThreshold]
    filteredFlux.drop(columns=['absFlux'], inplace=True)
    print(f'Number of reactions with absolute flux > {filterThreshold}: {len(filteredFlux)}')
    filteredQuantileVal = filteredFlux[f'{condition}Flux'].quantile(0.95)
    print(f'95% quantile of filtered fluxes: {filteredQuantileVal}')
    # a) heatmap of ordered filtered fluxes
    # order by flux
    filteredFlux = filteredFlux.sort_values(by=f'{condition}Flux', ascending=False)

    # set vmin vmax accordingly qanVal
    fig,ax = plt.subplots(figsize=(10, 20))
    sns.heatmap(filteredFlux[[f'{condition}Flux']], cmap='viridis', vmin=-5, vmax=5)
    # label y axis with "reactionIds"
    ax.set_yticklabels(filteredFlux['reactionId'], rotation=0)
    plt.show()
    fig.savefig(f'{outputDir}{condition}_filteredFlux_heatmap_{modelPrefix}.png', dpi=300)

def main():
    # read model
    model = read_sbml_model(_modelPath)

    # load flux frame
    fluxDf = pd.read_csv(_fluxDf)
    # heatmap of filtered flux
    filteredFlux(fluxDf, 0.3, _outputDir, _modelPrefix, _condition)
    # heatmap for summarised flux for each pathway

    # heatmap for pathways used

        # Pathway + Reaction heatmap
        # Heatmap for each pathway

        # heatmap for overview on pathways used (pathways x reactions)
    
    # correlation between number of reactions per group and average pathway flux
