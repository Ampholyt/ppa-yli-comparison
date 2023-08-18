# import 
import cobra
import numpy as np # load efm data
import matplotlib.pyplot as plt
import pandas as pd # upset plot
import sys # append path
from upsetplot import plot # upset plot
import os

sys.path.append('../scripts/')
import helperFunction as hf
# read config (for model path)
config = hf.load_config()

# read data
coefs = np.load(config['experiments']['efm_decomposition']['coefs'])
efms = np.load(config['experiments']['efm_decomposition']['efms'])
optimal_vector = np.load(config['experiments']['efm_decomposition']['optimal_vector'])

# set names of EFMs
_efm_names = [f'EFM {i}' for i in range(1, len(efms)+1)] 

# get number of active reactions in each EFM
_rxn_number = []
for efm in efms:
    _rxn_number.append(len(efm[np.where(efm!=0)])) 
# [336, 342, 340, 340, 2, 3, 3, 3]


if all(np.round(np.dot(efms.T,coefs),6) == np.round(optimal_vector,6)):
    print('Imported decomposition is equal to decomposition')


def scaled_efm(efm, coef):
    '''scale efms and their coefficients to the same range'''
    efm_scaled = efm/np.linalg.norm(efm)
    coef_scaled = coef * np.linalg.norm(efm)
    return efm_scaled, coef_scaled

def store_scaled_coefficients(outpath = '../results/EFM_decomposition/scaled_coefs.png'):
    '''store scaled coefficients plot'''
    scaled_efms = []
    scaled_coefs = []
    for i in range(8):
        efm_scaled, coef_scaled = scaled_efm(efms[i], coefs[i])
        scaled_efms.append(efm_scaled)
        scaled_coefs.append(coef_scaled)

    # plot scaled coefficients as bar plot with x-axis = EFMs annotate values
    plt.xlabel('EFMs')
    plt.ylabel('Scaled coefficients')
    for i in range(8):
        plt.bar(_efm_names[i], scaled_coefs[i])
        plt.annotate(round(scaled_coefs[i], 2), xy=(i, scaled_coefs[i]), ha='center', va='bottom')
    # store plot as png
    plt.savefig(outpath, dpi=300, bbox_inches='tight')

def plot_active_reactions(outpath='../results/EFM_decomposition/active_reactions.png'):
    '''plot number of active reactions in each EFM'''
    # inspired by # Source: https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 5))
    fig.subplots_adjust(hspace=0)  # adjust space between axes

    # same data on both axes
    ax1.bar(_efm_names, _rxn_number)
    ax2.bar(_efm_names, _rxn_number)

    # EFM 1-4
    first_end = 3.8
    second_start = 250
    second_end = 350

    # set value ranges
    ax1.set_ylim(second_start, second_end) # upper half
    ax2.set_ylim(0, first_end)  # lower half

    ax1.spines.bottom.set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop=False)
    ax2.spines.top.set_visible(False)
    ax2.xaxis.tick_bottom()

    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

    # # Add the values to the bars
    for i, v in enumerate(_rxn_number):
        ax1.text(i, v, str(v), ha='center', va='bottom')
        ax2.text(i, v, str(v), ha='center', va='bottom')

    # Add y-axis label for the whole figure
    fig.text(0.03, 0.5, 'Number of active reactions (non-zero fluxes)', va='center', rotation='vertical')
    plt.xlabel('EFMs')
    
    # store plot as png
    plt.savefig(outpath, dpi=300)

def get_rxns_from_indices(rxn_ids, indices):
    '''returns list of reactions with given indices'''
    rxns = []
    for i in indices:
        rxns.append(rxn_ids[i])
    return rxns

def upset_plot_first_four_efms(iLC915_model, outpath='../results/EFM_decomposition/upset_plot_EFM1-4.png'):
    '''plot upset plot for first four EFMs'''
    # combine EFMs and reaction IDs

    # get reaction IDs
    _rxm_ids = [rxn.id for rxn in iLC915_model.reactions]

    # get active reactions and indices from EFMs
    _efm_active_rxns = []
    _efm_active_index = []

    for efm in efms:
        non_zero_tuple = np.where(efm!=0)
        active_rxns = get_rxns_from_indices(_rxm_ids, non_zero_tuple[0])
        _efm_active_index.append(non_zero_tuple[0])
        _efm_active_rxns.append(active_rxns)


    rxn_set_1 = set(_efm_active_rxns[0])
    rxn_set_2 = set(_efm_active_rxns[1])
    rxn_set_3 = set(_efm_active_rxns[2])
    rxn_set_4 = set(_efm_active_rxns[3])


    set_names = ['EFM 1', 'EFM 2', 'EFM 3', 'EFM 4']
    all_elems = rxn_set_1.union(rxn_set_2).union(rxn_set_3).union(rxn_set_4)
    df = pd.DataFrame([[e in rxn_set_1, e in rxn_set_2, e in rxn_set_3, e in rxn_set_4] for e in all_elems], columns = set_names)
    df_up = df.groupby(set_names).size()
    plot(df_up, orientation='horizontal', sort_by='cardinality', show_counts=True)
    # store plot as png
    plt.savefig(outpath, dpi=300, bbox_inches='tight')


def main():
    '''main function'''
    print('EFM decomposition of iLC915 started ...')
    outpath_scaled = 'results/EFM_decomp_iLC915/scaled_coefs.png'
    outpath_active = '../results/EFM_decomposition/active_reactions.png'
    outpath_upset = '../results/EFM_decomposition/upset_plot_EFM1-4.png'
    
    if config['general']['snakemake']:
        outpath_scaled = snakemake.output[0]
        outpath_active = snakemake.output[1]
        outpath_upset = snakemake.output[2]
        
    # scaled coefficients
    store_scaled_coefficients(outpath_scaled)
    # active reactions
    plot_active_reactions(outpath_active)
    # upset plot
    iLC915_model = cobra.io.read_sbml_model(config['models']['ppaiLC915'])
    upset_plot_first_four_efms(iLC915_model, outpath_upset)
    print('EFM decomposition of iLC915 finished.')


if __name__ == '__main__':
    main()