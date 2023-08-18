# imports
import cobra
import yaml
from cobra.flux_analysis import pfba
# plotting
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


import sys # append path

sys.path.append('../scripts/')
import helperFunction as hf

config = hf.load_config()

_snakemake = config['general']['snakemake']
_show_plot = config['general']['show_plot']
_seperator = config['seperator']
_use_fructose = False
_store_results = config['general']['store_results']
_condition_name = config['general']['condition_name']
_table_out_path = f'../results/FBA_results/{_condition_name}_simulated_growth_rates.csv'
_outpath = f'../results/FBA_results/yli_glucose_growth_comparison.png'
if _snakemake: 
    _outpath = snakemake.output[0]
    _table_out_path = snakemake.output[1]
_add_maintenance = True
x_limit = 0.31 # default: 0.31
y_limit = 0.31 # default: 0.31
legend_location = 'upper left' # default: 'upper left' other: 'best'
# _NGAM = (1,1) # like in iMK735
# _NGAM = (5.03, 5.03) # like in literature iYLI647
_NGAM = (7.8625, 7.8625) # like in literature iYli21 + iYli_2.0


# experimental data:
# all 
yli_glucose_uptake = [2.43, 0.64, 0.95, 0.33, 0.72, 0.98, 2.09] # from 2 different resources (xu)
yli_experimental_growth = [0.26, 0.048, 0.08, 0.03, 0.07, 0.1, 0.2]

# # only W29
# yli_glucose_uptake = [2.43, 0.64, 0.61] # from 2 different resources (Xu, Guo, Lazar)
# yli_experimental_growth = [0.26, 0.048, 0.047]

# everthing is normoxic for yli (27.03.2023)
# # experimental data from Xu et al. 2020 Comparison: https://link.springer.com/article/10.1007/s12257-019-0208-1#Sec1 supplementary table 1
# ref; Uptake Glucose; Growth rate
# [1]     0.95             0.08
# [2]     0.33             0.03
# [2]     0.72             0.07
# [2]     0.98             0.1
# [2]     2.09             0.2
# [4]     2.46             0.24 (suspicious because no information about oxygen concentration provided)
## yli: https://www.sciencedirect.com/science/article/pii/S200103702200174X?via%3Dihub



# load models 
# iYli21:
print('iYli21_model')
iYli21_model_path = config['models']['yli21']
iYli21_model = cobra.io.read_sbml_model(iYli21_model_path)
# set objective
iYli21_model.objective = 'biomass_C'

glucose_exchange_reaction = 'R1070'
fructose_exchange_reaction = 'R1065'
biomass_reaction = 'biomass_C'
maintenance_reaction = 'xMAINTENANCE'


# Set glucose reaction to -2.43 mmol/h, 0.64 and 0.61 (based on Guo et al. 2022)
iYli21_simulation_growth = []
for condition in yli_glucose_uptake:
    with iYli21_model as model:
        if _use_fructose:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (0, 1000)
            dGlu = model.reactions.get_by_id(fructose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        else:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        if _add_maintenance:
            maintenacne_reaction = model.reactions.get_by_id(maintenance_reaction)
            maintenacne_reaction.bounds = _NGAM
        try:
            pfba_solution = pfba(model)
            growthRate = pfba_solution.fluxes[biomass_reaction]
        except:
            growthRate = 0
        iYli21_simulation_growth.append(growthRate)
        print(f'Growth rate on {condition} mmol/h {_condition_name}: {growthRate}')


# iYali4:
print('iYali4_model')
iYali4_model = cobra.io.read_sbml_model(config['models']['yli4_corr'])

# set objective
iYali4_model.objective = 'biomass_C'

# growth on glucose
glucose_exchange_reaction = '1714'
fructose_exchange_reaction = '1709'
biomass_reaction = 'biomass_C'

# maintenance
maintenance_reaction = 'xMAINTENANCE'

# Set glucose reaction to -2.43 mmol/h, 0.64 and 0.61 (based on Guo et al. 2022)
iYali4_simulation_growth = []
for condition in yli_glucose_uptake:
    with iYali4_model as model:
        if _use_fructose:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (0, 1000)
            dGlu = model.reactions.get_by_id(fructose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        else:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        if _add_maintenance:
            maintenacne_reaction = model.reactions.get_by_id(maintenance_reaction)
            maintenacne_reaction.bounds = _NGAM
        try:
            pfba_solution = pfba(model)
            growthRate = pfba_solution.fluxes[biomass_reaction]
        except:
            growthRate = 0
        iYali4_simulation_growth.append(growthRate)
        print(f'Growth rate on {condition} mmol/h {_condition_name}: {growthRate}')


# iMK735:
print('iMK735_model')
iMK735_model_path = config['models']['yliMK735_corr']
iMK735_model = cobra.io.read_sbml_model(iMK735_model_path)


# set objective
iMK735_model.objective = 'biomass_C'

glucose_exchange_reaction = 'EX_glc(e)'
fructose_exchange_reaction = 'EX_fru(e)'
biomass_reaction = 'biomass_C'

# maintenance
maintenance_reaction = 'ATPM'

iMK735_simulation_growth = []
for condition in yli_glucose_uptake:
    with iMK735_model as model:
        if _use_fructose:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (0, 1000)
            dGlu = model.reactions.get_by_id(fructose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        else:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        if _add_maintenance:
            maintenacne_reaction = model.reactions.get_by_id(maintenance_reaction)
            maintenacne_reaction.bounds = _NGAM
        try:
            pfba_solution = pfba(model)
            growthRate = pfba_solution.fluxes[biomass_reaction]
        except:
            growthRate = 0
        iMK735_simulation_growth.append(growthRate)
        print(f'Growth rate on {condition} mmol/h {_condition_name}: {growthRate}')


# iNL895:
print('iNL895_model')
iNL895_model_path = config['models']['yliNL895_corr']
iNL895_model = cobra.io.read_sbml_model(iNL895_model_path)

# set objective
iNL895_model.objective = 'biomass_C'

# maintenance: not existing/ not found


# Set glucose reaction to -2.43 mmol/h, 0.64 and 0.61 (based on Guo et al. 2022)
biomass_reaction = 'biomass_C' # 'biomass_C' # 'r_1814' # 'r_021_xxx'
glucose_exchange_reaction = 'r_51_exchange'
iNL895_simulation_growth = []
for condition in yli_glucose_uptake:
    with iNL895_model:
        dGlu = iNL895_model.reactions.get_by_id(glucose_exchange_reaction)
        dGlu.bounds = (-1000, condition)
        try:
            pfba_solution = pfba(iNL895_model)
            growthRate = pfba_solution.fluxes[biomass_reaction]
        except:
            growthRate = 0
        iNL895_simulation_growth.append(growthRate)
        print(f'Growth rate on {condition} mmol/h {_condition_name}: {growthRate}')

# iYli_2.0: 
print('iyli_2.0_model')
iyli2_0_model_path = config['models']['yli2.0_corr']
iyli2_0_model = cobra.io.read_sbml_model(iyli2_0_model_path)


# set objective
iyli2_0_model.objective = 'biomass_C'

# growth on glucose
glucose_exchange_reaction = 'R1294'
fructose_exchange_reaction = 'R1261'
biomass_reaction = 'biomass_C'

# maintenance
maintenance_reaction = 'R0542'

# Set glucose reaction to -2.43 mmol/h, 0.64 and 0.61 (based on Guo et al. 2022)
iyli2_0_simulation_growth = []
for condition in yli_glucose_uptake:
    with iyli2_0_model as model:
        if _use_fructose:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (0, 1000)
            dGlu = model.reactions.get_by_id(fructose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        else:
            dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
            dGlu.bounds = (-condition, 1000)
        if _add_maintenance:
            maintenacne_reaction = model.reactions.get_by_id(maintenance_reaction)
            maintenacne_reaction.bounds = _NGAM
        try:
            pfba_solution = pfba(model)
            growthRate = pfba_solution.fluxes[biomass_reaction]
        except:
            growthRate = 0
        iyli2_0_simulation_growth.append(growthRate)
        print(f'Growth rate on {condition} mmol/h {_condition_name}: {growthRate}')

# iYLI647:
print('iYLI647_model')
iyli647_model_path = config['models']['yli647_corr']
iyli647_model = cobra.io.read_sbml_model(iyli647_model_path)

# set objective
iyli647_model.objective = 'biomass_C'

# Generate data for the plot
glucose_exchange_reaction = 'EX_glc(e)'
biomass_reaction = 'biomass_C'

# maintenance
maintenance_reaction = 'ATPM'

# Set glucose reaction to -2.43 mmol/h, 0.64 and 0.61 (based on Guo et al. 2022)
iyli647_simulation_growth = []
for condition in yli_glucose_uptake:
    with iyli647_model as model:
        dGlu = model.reactions.get_by_id(glucose_exchange_reaction)
        dGlu.bounds = (-condition, 1000)
        if _add_maintenance:
            maintenacne_reaction = model.reactions.get_by_id(maintenance_reaction)
            maintenacne_reaction.bounds = _NGAM
        try:
            pfba_solution = pfba(model)
            growthRate = pfba_solution.fluxes[biomass_reaction]
        except:
            growthRate = 0
        iyli647_simulation_growth.append(growthRate)
        print(f'Growth rate on {condition} mmol/h {_condition_name}: {growthRate}')

# dataframe of experimental and simulated growth rates
# generate dataframe with experimental results
simulated_growth_df = pd.DataFrame({
    f'{_condition_name}_uptake': yli_glucose_uptake,
    'experimental_growth': yli_experimental_growth,
    'iYli21_model': iYli21_simulation_growth,
    'iYali4(cor)_model': iYali4_simulation_growth,
    'iMK735(cor)_model': iMK735_simulation_growth,
    'iNL895(cor)_model': iNL895_simulation_growth,
    'iYli_2_0(cor)_model': iyli2_0_simulation_growth,
    'iYLI647(cor)_model': iyli647_simulation_growth
    })

# write dataframe to csv file
if _store_results:        
    simulated_growth_df.to_csv(_table_out_path, sep=_seperator, index=False)
    print(f'Wrote simulated growth rates to {_table_out_path}')

# plot experimental and simulated growth rates
print('\nPlotting experimental and simulated growth rates ...\n')

# plot
fig, ax = plt.subplots(figsize=(5,5))
sns.scatterplot(x=yli_experimental_growth, y=iYli21_simulation_growth, ax=ax, color='blue', marker="o", label='iYli21')
sns.scatterplot(x=yli_experimental_growth, y=iYali4_simulation_growth, ax=ax, color='orange', marker="o", label='iYali4(cor)')
sns.scatterplot(x=yli_experimental_growth, y=iNL895_simulation_growth, ax=ax, color='purple', marker="o", label='iNL895(cor)')
sns.scatterplot(x=yli_experimental_growth, y=iyli2_0_simulation_growth, ax=ax, color='black', marker="o", label='iYli2.0(cor)')
sns.scatterplot(x=yli_experimental_growth, y=iyli647_simulation_growth, ax=ax, color='brown', marker="o", label='iYLI647(cor)')
sns.scatterplot(x=yli_experimental_growth, y=iMK735_simulation_growth, ax=ax, color='yellow', marker="o", label='iMK735(cor)')
ax.set_xlabel('Experimental growth rate [1/h]')
ax.set_ylabel('Simulated growth rate [1/h]')
ax.set_title('All corrected models from Chen et al. 2016 and iYli21 on glucose\niYli21 vs iYali4 vs iNL895 vs iYli2.0 vs iYLI647 vs iMK735 vs experimental data')

if _add_maintenance:
    ax.set_title('All corrected models from Chen et al. 2016 and iYli21 on glucose\niYli21 vs iYali4 vs iNL895 vs iYli2.0 vs iYLI647 vs iMK735 vs experimental data' + '\n' + f'with maintenance: {_NGAM[0]} ATP mmol/gDW/h')

ax.legend(loc=legend_location)

# limit x and y axis
plt.ylim([0, y_limit])
plt.xlim([0, x_limit])

# plot diagonal line
x = np.linspace(0, x_limit, 100)
ax.plot(x, x, color='red', linestyle='--')

# save plot
# make sure the output path exists
os.makedirs(os.path.dirname(_outpath), exist_ok=True)
plt.savefig(_outpath)

if _show_plot:
    plt.show()

########################  Plotting as in XU et al. 2020  #########################################################################

# # only xu experimental data
# yli_glucose_uptake = [0.95, 0.33, 0.72, 0.98, 2.09, 0.61, 0.64, 2.46]
# yli_experimental_growth = [0.08, 0.03, 0.07, 0.1, 0.2, 0.047, 0.048, 0.24]

# sns.scatterplot(x=yli_experimental_growth, y=iyli647_simulation_growth, ax=ax, color='magenta', marker="<", label='iYLI647')
# sns.scatterplot(x=yli_experimental_growth, y=iYali4_simulation_growth, ax=ax, color='black', marker="s", label='iYali4')
# sns.scatterplot(x=yli_experimental_growth, y=iyli2_0_simulation_growth, ax=ax, color='green', marker="v", label='iYli2.0')
# sns.scatterplot(x=yli_experimental_growth, y=iMK735_simulation_growth, ax=ax, color='red', marker="o", label='iMK735')
# sns.scatterplot(x=yli_experimental_growth, y=iNL895_simulation_growth, ax=ax, color='blue', marker="^", label='iNL895')

