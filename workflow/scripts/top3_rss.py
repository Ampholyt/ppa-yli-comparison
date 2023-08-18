# script:

# import 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt #plot 
import os # path.join
import scipy.stats as stats # z score normalization

import itertools # plot of the top 3 perfoming models

# config preperation
import sys # append path

sys.path.append('../scripts/')
import helperFunction as hf

config = hf.load_config()

def plot_RSS(value_dict, conditions, title, y_label, color_list = []):
    '''Plot the RSS for each model given by the value_dict
    e.g. value_dict = {'iYali4_model': [0.0036525137174046765],...}
    @params: min and max values are used for limiting the y-axis
    @trick: (is extendable) is able to plot conditions for different models'''
    min_value = min(value_dict.values())[0]
    max_value = max(value_dict.values())[-1]
    
    x = np.arange(len(conditions))  # the label locations
    width = 0.01  # the width of the bars
    multiplier = 0
    fig, ax = plt.subplots(layout='constrained')

    
    idx = 0
    for attribute, measurement in value_dict.items():
        offset = width * multiplier
        print(attribute.split('_')[0])
        if color_list != []:
            rects = ax.bar(x + offset, measurement, width, label=attribute.split('_')[0], color=color_list[idx])
        else:
            rects = ax.bar(x + offset, measurement, width, label=attribute.split('_')[0])
        # ax.bar_label(rects, padding=3)
        multiplier += 1
        idx += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.set_xticks(x + width, conditions)
    ax.legend(loc='upper left', ncol=3)
    ax.set_ylim(min_value - 0.1 * abs(min_value), max_value + 0.1 * abs(max_value))
    return fig

def saveFig(outpath, figure):
    """Save figure in outpath"""
    figure.savefig(outpath, dpi=300, bbox_inches='tight')
    if config['general']['show_plot']:
        plt.show()

def main():
    # load dataframe with experimental results
    condition_name = config['general']['condition_name']
    growth_rates = pd.read_csv(f'../results/FBA_results/{condition_name}_simulated_growth_rates.csv', sep=config['seperator'])

    # for each column except experimental growth rate and uptake compute residual sum of squares
    experimental_col = [col for col in growth_rates.columns if "experimental" in col][0]
    simulation_cols = [col for col in growth_rates.columns if "experimental" not in col and "uptake" not in col]

    # get dict with col_name: rss
    rss_dict = {}

    for col_name in simulation_cols:
        rss_dict[col_name] = np.sum(np.square(growth_rates[col_name] - growth_rates[experimental_col]))
        # print(f'{col_name}: \n\tresidual sum of squares is : {rss_dict[col_name]} ')

    # sort dict by rss
    rss_dict = {k: [v] for k, v in sorted(rss_dict.items(), key=lambda item: item[1])}

    top_3_rss = dict(itertools.islice(rss_dict.items(), 3))
    conditions = ["Glucose Experiments\n(from Chan et al. 2017 and Guo et al. 2022)"]

    ylabel = 'Residual sum of squares (RSS)'
    title = 'Model quality using Residual sum of squares'
    fig = plot_RSS(top_3_rss, conditions, title, ylabel, color_list=['red', 'orange', 'green'])

    # Save figure
    _snakemake = config['general']['snakemake']
    outpath = os.path.join(config['experiments']['yli_growth'], 'yli_top3_rss.png')
    if _snakemake:
        outpath = snakemake.output[0]
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    if config['general']['show_plot']:
        plt.show()

    saveFigure(outpath, fig)


    # code from stackoverflow (https://stackoverflow.com/questions/62104874/multiple-bars-in-one-bar-plot)
    # plot of all models

    _use_normalization = False

    conditions = ["Glucose Experiments"]

    if _use_normalization:
        rss_df = pd.DataFrame.from_dict({'model': rss_dict.keys(), 'rss_raw': [val[0] for val in rss_dict.values()]})
        rss_df['rss_norm'] = stats.zscore(rss_df['rss_raw'])
        zipped_dict = dict(zip(rss_df['model'], rss_df['rss_norm']))
        rss_dict = {key: [value] for key, value in zipped_dict.items()}

    ylabel = 'Residual sum of squares (RSS)'
    title = 'Model quality using Residual sum of squares'
    fig = plot_RSS(rss_dict, conditions, title, ylabel, color_list=['red', 'orange', 'green', 'blue', 'purple', 'brown'])

    if _snakemake:
        outpath = snakemake.output[1]


if __name__ == '__main__':
    main()

