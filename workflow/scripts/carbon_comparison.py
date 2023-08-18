## Skript which uses the information of the different usable carbon sources and their comparison between the two species Ppa and Yli (FBA)
# imports
import cobra
from cobra.flux_analysis import pfba # parsimonious FBA
import matplotlib.pyplot as plt
import numpy as np
import sys # append path

sys.path.append('../scripts/')
import helperFunction as hf

config = hf.load_config()

# parameters
_use_precomputed = False
_verbose = config['general']['verbose']
_snakemake = config['general']['snakemake']
_show_plot = config['general']['show_plot']
_experiment_number = 8
_carbon_sources = ['Alanine (L-Alanine)', 'Methanol', 'Oleic acid', 'Glucose', 'Fructose', 'Trehalose', 'Sorbitol', 'Glycerol']
_number_of_carbons = [3, 1, 18, 6, 6, 12, 6, 3]
_ppa_maintenance = 2.81 # according to iMT1026v3
_add_maintenance = True
_yli_maintenance = 7.8625

# load models (iYli21, iYali4 (corrected), iMT1026v3 (read_write), iLC915 (MODEL1507))
iYli21 = cobra.io.read_sbml_model(config['models']['yli21'])
iYali4 = cobra.io.read_sbml_model(config['models']['yli4_corr'])
iMT1026v3 = cobra.io.read_sbml_model(config['models']['ppa1026v3'])
iLC915 = cobra.io.read_sbml_model(config['models']['ppaiLC915'])

# # Alanine, Methanol, Oleic acid, Glucose, Fructose, Trehalose, Sorbitol, Glycerol, CO2 productio, O2 production

yli21_carbon = ['R1196', '-', 'R1440', 'R1070', 'R1065', 'R1013', 'R1209', 'R1141'] # , 'R1030', 'R1287'
yali4_carbon = ['1873', '-', '2189', '1714', '1709', '1650', '1712', '1808'] # , '1672', '1992'
iLC915_carbon = ['r1126', 'r1158', 'r1161', 'r1145', 'r1144', 'r1173', 'r1170', 'r1148', 'r1137', 'r1160']
iMT1026v3_carbon = ['Ex_ala_L', 'Ex_meoh', '-', 'Ex_glc_D', 'Ex_fru', 'Ex_tre', 'Ex_sbt_D', 'Ex_glyc'] # , 'Ex_co2', 'Ex_o2'

def prepare_iLC915():
    '''Prepare iLC915 model requried because internal cycles'''
    if _verbose:
        print('iLC915 model: ')
    # prepare iLC915: needed because of internal cycles and active reactions
    # lower and upper bounds of several reactions should be set to zero according to supplementary material of 
    constraints_2_zero = ['r66','r910','r1104','r239','r111','r106','r490','r791','r243','r252','r253','r307','r308','r404','r405','r1320','r639','r640','r641','r642','r649','r650','r651','r652','r645','r646','r643','r644','r653','r654','r655','r656','r534'] 
    for rxn_id in constraints_2_zero:
        rxn = iLC915.reactions.get_by_id(rxn_id)
        rxn.bounds = (0,0)
    # constraint formulate uptake because of not meaningful results
    formulate_uptake = iLC915.reactions.get_by_id('r1143')
    formulate_uptake.bounds = (0,0)
    # important carbon sources (biotin, CO2, urea)
    important_carbon = ['r1132','r1137','r1177']

    # all carbon containing exchange reactions
    carbon_exchange_rxns = ['r1122', 'r1123', 'r1124', 'r1126', 'r1127', 'r1129', 'r1130', 'r1131', 'r1132', 'r1134', 'r1135', 'r1137', 'r1138', 'r1139', 'r1140', 'r1141', 'r1144', 'r1145', 'r1146', 'r1147', 'r1148', 'r1149', 'r1151', 'r1152', 'r1153', 'r1154', 'r1155', 'r1156', 'r1157', 'r1158', 'r1161', 'r1162', 'r1163', 'r1165', 'r1167', 'r1168', 'r1170', 'r1172', 'r1173', 'r1174', 'r1175', 'r1176', 'r1177', 'r1178']

    # set all carbon exchange reactions to 0 except important ones
    for rxn_id in carbon_exchange_rxns:
        if rxn_id not in important_carbon:
            rxn = iLC915.reactions.get_by_id(rxn_id)
            rxn.bounds = (0,0)

    # check definition of important rxns 
    important_rxns = ['r1133', 'r1141', 'r1145', 'r1144', 'r1148', 'r1160', 'r1137', 'r1188']

    if _verbose: 
        for rxn_id in important_rxns:
            rxn = iLC915.reactions.get_by_id(rxn_id)
            print(rxn.reactants, rxn.products, rxn.bounds)
            print(rxn_id, rxn.name, hf.formulaWithNames(rxn), 'rxn coefficient: ', rxn.get_coefficient(list(rxn.metabolites.keys())[0]))

    biomass_rxn = 'r1339' # growth
    # set objective: 
    iLC915.objective = biomass_rxn

    # set maintenance reaction
    maintenance = 'r1188'
    maintenance_rxn = iLC915.reactions.get_by_id(maintenance)
    maintenance_rxn.bounds = (1.2434,1.2434) # according to iLC915 (1.2434 * 2.26 = 2.81)
        
    try:
        iLC915.summary()
        raise ValueError('Model behaves not as expected')

    except: 
        if _verbose:
            print('model works as expected')

def revert_ilc915():
    iLC915_carbon = ['r1126', 'r1158', 'r1161', 'r1145', 'r1144', 'r1173', 'r1170', 'r1148'] # , 'R1030', 'R1287'
    for rxn_id in iLC915_carbon:
        if rxn_id == '-':
            continue
        rxn = iLC915.reactions.get_by_id(rxn_id)
        rxn.bounds = (0,0)

def compute_iLC915():
    """Compute growth data for Yli model iLC915"""
    biomass_rxn = 'r1339' # growth
    iLC915_growth = []
    iLC915_biomass_ylied = []
    if not _use_precomputed:
        for i in range(len(_number_of_carbons)):
            uptake_carbon = 1000/_number_of_carbons[i]
            carbon_rxn = iLC915.reactions.get_by_id(iLC915_carbon[i])
            carbon_rxn.bounds = (-uptake_carbon, uptake_carbon)
            print(carbon_rxn.name)
            try: 
                sol = pfba(iLC915)
            except:
                print('infeasible')
                iLC915_growth.append(0)
                iLC915_biomass_ylied.append(0)
                revert_ilc915()
                continue
            growth_rate = round(sol.fluxes[biomass_rxn],4)
            print(growth_rate) 
            biomass_yield = round(growth_rate/uptake_carbon,4)
            print(biomass_yield)
            iLC915_growth.append(growth_rate)
            iLC915_biomass_ylied.append(biomass_yield)
            revert_ilc915()
    else:
        iLC915_growth = [13.1491, 8.8747, 18.5102, 16.6677, 0, 16.67, 17.859, 18.8556]
        iLC915_biomass_ylied = [0.0394, 0.0089, 0.3332, 0.1, 0, 0.2, 0.1072, 0.0566]
    return iLC915_growth, iLC915_biomass_ylied

def prepare_iMT1026v3():
    """Prepare iMT1026v3: no growth possible after preparation"""
    # prepare the model: no growth possible: missing carbon source
    # show interesting reactions
    model_name = 'iMT1026v3'
    maintenance = 'ATPM'
    biomass_reaction = 'growth' # default: 'growth' biomass
    glu_ex_rxn = 'Ex_glc_D' # glucose
    glycerol_ex_rxn = 'Ex_glyc' # glycerol
    o2_ex_rxn = 'Ex_o2' # O2
    co2_ex_rxn = 'Ex_co2' # CO2
    fructose_ex_rxn = 'Ex_fru' # fructose
    ethanol_ex_rxn = 'Ex_etoh' # ethanol
    arabitol_ex_rxn = 'Ex_abt_D' # Arabitol/arabinitol

    # check definition of interesting rxns
    interesting_rxns = ['growth', 'Ex_glc_D', 'ATPM', 'Ex_glyc', 'Ex_o2', 'Ex_co2', 'Ex_fru', 'Ex_etoh', 'Ex_abt_D']

    if _verbose:
        for rxn_id in interesting_rxns:
            rxn = iMT1026v3.reactions.get_by_id(rxn_id)
            print(rxn.reactants, rxn.products, rxn.bounds)
            print(rxn_id, rxn.name, hf.formulaWithNames(rxn), 'rxn coefficient: ', rxn.get_coefficient(list(rxn.metabolites.keys())[0]), rxn.bounds)

    # remove glyc reaction as carbon source
    glyc_rxn = iMT1026v3.reactions.get_by_id('Ex_glyc')
    glyc_rxn.bounds = (0, 1000)

    # add o2 bounds
    o2_bounds = (-1000, 1000)
    Ex_o2 = iMT1026v3.reactions.get_by_id(o2_ex_rxn)
    Ex_o2.bounds = o2_bounds

    # set objective
    iMT1026v3.objective = biomass_reaction

    # set maintenance 
    maintenance_rxn = model.reactions.get_by_id(maintenance)
    maintenance_rxn.bounds = (_ppa_maintenance,_ppa_maintenance)
 
    try:
        iMT1026v3.summary() # not feasible => all good
        # raise error
        raise ValueError('Model behaves not as expected')
    except:
        print('Model behaves as expected')

def compute_iMT1026v3():
    """Compute growth data for Ppa model iMT1026v3"""
    if _verbose:
        print('iMT1026v3_model')
    if _use_precomputed:
        iMT1026v3_simulation_growth = [0, 9.226, 10.9062, 15.0914, 15.0914, 15.0952, 15.0914, 14.8574]
        iMT1026v3_simulation_biomass_yield = [0, 0.0092, 0.1963, 0.0905, 0.0905, 0.1811, 0.0905, 0.0446]
        return iMT1026v3_simulation_growth, iMT1026v3_simulation_biomass_yield

    maintenance = 'ATPM'
    biomass_reaction = 'growth' # default: 'growth' biomass
    # set objective
    iMT1026v3.objective = biomass_reaction

    # growth results
    iMT1026v3_simulation_growth = []
    iMT1026v3_simulation_biomass_yield = []
    iMT1026v3_simulation_co2_prod = []
    iMT1026v3_simulation_o2_prod = []
    
    for exp_idx in range(_experiment_number):
        if iMT1026v3_carbon[exp_idx] == '-':
            print('No carbon source')
            iMT1026v3_simulation_growth.append(0)
            iMT1026v3_simulation_biomass_yield.append(0)
            iMT1026v3_simulation_co2_prod.append(0)
            iMT1026v3_simulation_o2_prod.append(0)
            continue
        with iMT1026v3 as model:
            # set carbon source accroding to uptake reaction and set bounds according to carbon number
            carbon_rxn = model.reactions.get_by_id(iMT1026v3_carbon[exp_idx]) 
            uptake_rate = 1000/_number_of_carbons[exp_idx]
            carbon_rxn.bounds = (-uptake_rate, uptake_rate) # rxn is defined as: "glc_D_e -->"
            print(f'{iMT1026v3_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h')

            # set maintenance
            maintenance_rxn = model.reactions.get_by_id(maintenance)
            maintenance_rxn.bounds = (_ppa_maintenance,_ppa_maintenance)

            # get solution
            try: 
                solution = pfba(model)
                growth = round(solution.fluxes[biomass_reaction],4)
                biomass_yield = round(growth/uptake_rate, 4)
                co2_prod = round(solution.fluxes[co2_ex_rxn],4)
                o2_prod = round(solution.fluxes[o2_ex_rxn],4)
            except:
                print('growth not possible')
                growth = 0
                biomass_yield = 0
                co2_prod = 0
                o2_prod = 0
            print(f'Growth rate on {iMT1026v3_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h: {growth}')
            print(f'Biomass yield on {iMT1026v3_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h: {biomass_yield}')
            print(f'CO2 production on {iMT1026v3_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h: {co2_prod}')
            print(f'O2 production on {iMT1026v3_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h: {o2_prod}')

            # store growth
            iMT1026v3_simulation_growth.append(growth)
            iMT1026v3_simulation_biomass_yield.append(biomass_yield)
            # store co2 production
            iMT1026v3_simulation_co2_prod.append(co2_prod)
            # store o2 production
            iMT1026v3_simulation_o2_prod.append(o2_prod)
    return iMT1026v3_simulation_growth, iMT1026v3_simulation_biomass_yield


def revert_yli21():
    yli21_carbon = ['R1196', '-', 'R1440', 'R1070', 'R1065', 'R1013', 'R1209', 'R1141'] # , 'R1030', 'R1287'
    for rxn_id in yli21_carbon:
        if rxn_id == '-':
            continue
        rxn = iYli21.reactions.get_by_id(rxn_id)
        rxn.bounds = (0,0)

def compute_yli21():
    if _verbose:
        print('iYli21_model')

    if _use_precomputed:
        iYli21_simulation_growth = [18.4807, 0, 18.4804, 22.9298, 22.9298, 0, 0, 25.6036]
        iYli21_simulation_biomass_yield = [0.0554, 0.3326, 0.1376, 0.1376, 0, 0, 0.0768]
        return iYli21_simulation_growth, iYli21_simulation_biomass_yield
    # set objective
    iYli21.objective = 'biomass_C'
    glucose_exchange_reaction = 'R1070'
    fructose_exchange_reaction = 'R1065'
    biomass_reaction = 'biomass_C'
    maintenance_reaction = 'xMAINTENANCE'
    co2_ex_rxn = 'R1030'
    o2_ex_rxn = 'R1287'

    # remove glucose as carbon source
    revert_yli21()

    iYli21_simulation_growth = []
    iYli21_simulation_biomass_yield = []
    iYli21_simulation_o2 = []
    iYli21_simulation_co2 = []
    for exp_idx in range(_experiment_number):
        with iYli21 as model:
            try:
                model.summary()
                print('Model does not work as expected')
            except:
                print('Model works as expected')
                print('Simulation started...')
            if yli21_carbon[exp_idx] == '-':
                print('No carbon source')
                iYli21_simulation_growth.append(0)
                iYli21_simulation_co2.append(0)
                iYli21_simulation_o2.append(0)
                continue
            # set carbon source accroding to uptake reaction and set bounds according to carbon number
            carbon_rxn = model.reactions.get_by_id(yli21_carbon[exp_idx]) 
            uptake_rate = 1000/_number_of_carbons[exp_idx]
            carbon_rxn.bounds = (-uptake_rate, uptake_rate) # rxn is defined as: "glc_D_e -->"
            print(f'{yli21_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h')
            # set carbon source accroding to uptake reaction and set bounds according to carbon number

            if _add_maintenance:
                maintenacne_reaction = model.reactions.get_by_id(maintenance_reaction) # m141[c] --> m143[c] + m35[c] => lower bound > 0
                maintenacne_reaction.bounds = (_yli_maintenance, _yli_maintenance)
            try: 
                solution = pfba(model)
                growth = round(solution.fluxes[biomass_reaction],4)
                biomass_yield = round(growth/uptake_rate, 4)
                co2_prod = round(solution.fluxes[co2_ex_rxn],4)
                o2_prod = round(solution.fluxes[o2_ex_rxn],4)
            except:
                print('growth not possible')
                biomass_yield = 0
                growth = 0
                co2_prod = 0
                o2_prod = 0
            print(f'Growth rate on {yli21_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h: {biomass_yield}')
            print(f'CO2 production on {yli21_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h: {co2_prod}')
            print(f'O2 production on {yli21_carbon[exp_idx]} with {uptake_rate} mmol/gDW/h: {o2_prod}')
            iYli21_simulation_growth.append(growth)
            iYli21_simulation_biomass_yield.append(biomass_yield)
            iYli21_simulation_co2.append(co2_prod)
            iYli21_simulation_o2.append(o2_prod)
    return iYli21_simulation_growth, iYli21_simulation_biomass_yield

def revert_yali4():
    iYali4rxns = ['1873', '-', '2189', '1714', '1709', '1650', '1712', '1808']
    for id in iYali4rxns:
        if id == '-':
            continue
        rxn = iYali4.reactions.get_by_id(id)
        rxn.bounds = (0,0)

def compute_yali4():
    """Compute growth data for Yli model iYali4"""
    
    if _verbose:
        print('iYali4_model')
    if _use_precomputed:
        yali4_growth_rate = [18.9218, 0, 21.2723, 22.9068, 22.9068, 22.9097, 24.7162, 25.5554]
        yali4_biomass_yield = [0.0568, 0, 0.3829, 0.1374, 0.1374, 0.2749, 0.1483, 0.0767]
        return yali4_growth_rate, yali4_biomass_yield
    
    # set maintenance
    xMAINTENANCE = iYali4.reactions.get_by_id('xMAINTENANCE')
    xMAINTENANCE.bounds = (_yli_maintenance, _yli_maintenance)

    biomass_rxn = 'biomass_C'
    iYali4.objective = biomass_rxn

    yali4_biomass_yield = []
    yali4_growth_rate = []
    for exp_idx in range(_experiment_number):
        revert_yali4()
        uptake_carbon = 1000/_number_of_carbons[exp_idx]
        carbon_rxn_id = yali4_carbon[exp_idx]
        if carbon_rxn_id == '-':
            print('No carbon source')
            yali4_biomass_yield.append(0)
            yali4_growth_rate.append(0)
            continue
        # set glucose 
        carbon_rxn = iYali4.reactions.get_by_id(carbon_rxn_id)
        carbon_rxn.bounds = (-uptake_carbon, uptake_carbon)
        print(carbon_rxn.name)
        sol = pfba(iYali4)
        growth_rate = round(sol.fluxes[biomass_rxn],4)
        print(growth_rate) 
        yali4_growth_rate.append(growth_rate)
        biomass_yield = round(growth_rate/uptake_carbon,4)
        print(biomass_yield)
        yali4_biomass_yield.append(biomass_yield)
    return yali4_growth_rate, yali4_biomass_yield

def plot_carbon_sources():
    """Plot carbon sources using growth rates"""
    biomass_yield_yli21 = iYli21_biomass_ylied # [0.0554, 0.3326, 0.1376, 0.1376, 0, 0, 0.0768]
    growth_yli21 = iYli21_simulation_growth # [18.4807, 0, 18.4804, 22.9298, 22.9298, 0, 0, 25.6036]


    biomass_yield_yali4 = iYali4_biomass_yield
    growth_yali4 = iYali4_growth_rate
    # biomass_yield_yali4 = [0.0568, 0, 0.3829, 0.1374, 0.1374, 0.2749, 0.1483, 0.0767]
    # growth_yali4 = [18.9218, 0, 21.2723, 22.9068, 22.9068, 22.9097, 24.7162, 25.5554]

    growth_iLC915 = iLC915_growth # [13.1491, 8.8747, 18.5102, 16.6677, 0, 16.67, 17.859, 18.8556]
    biomass_yield_iLC915 = iLC915_biomass_ylied # [0.0394, 0.0089, 0.3332, 0.1, 0, 0.2, 0.1072, 0.0566]

    growth_iMT1026v3 = iMT1026v3_simulation_growth # [0, 9.226, 10.9062, 15.0914, 15.0914, 15.0952, 15.0914, 14.8574]
    biomass_yield_iMT1026v3 = iMT1026v3_simulation_biomass_yield # [0, 0.0092, 0.1963, 0.0905, 0.0905, 0.1811, 0.0905, 0.0446]

    # Set the width of the bars
    bar_width = 0.2

    # Calculate the position of the bars on the y-axis
    y_pos = np.arange(len(_carbon_sources))

    model_names = ['Y. lipolytica (iYali4 (cor))', 'Y. lipolytica (iYli21)', 'P. pastoris (iMT1026v3)', 'P. pastoris (iLC915)']

    # Set the figure size
    plt.figure(figsize=(8, 6))  # Adjust the values (width, height) as desired

    # Create the first subplot for growth of ppa (iLC915 vs iMT1026v3)
    plt.barh(y_pos + bar_width/0.9, growth_iMT1026v3, height=bar_width/1.5, label=model_names[2], color='turquoise')
    plt.barh(y_pos + bar_width/2.5, growth_iLC915, height=bar_width/1.5, label=model_names[3], color='blue')
    plt.barh(y_pos - bar_width/2.5, growth_yli21, height=bar_width/1.5, label=model_names[1], color='orange')
    plt.barh(y_pos - bar_width/0.9, growth_yali4, height=bar_width/1.5, label=model_names[0], color='red')

    # Add horizontal lines
    for i in range(len(_carbon_sources) - 1):
        plt.axhline(y=i+0.5, color='gray', linestyle='--', linewidth=1)

    # fix x axis limit to 30
    plt.xlim(0, 30)

    plt.xlabel('Growth rate (1/h)')
    plt.ylabel('Carbon Sources (1 Carbon/mol)')
    plt.title('')
    plt.yticks(y_pos, _carbon_sources)
    plt.legend(loc='lower right')
    fig1 = plt.gcf()
    if _show_plot:
        plt.show()
    # store the figure in good quality in ../results/FBA_results/growth_rate_plots/ folder with the name: Ppa_Yli_biomass_yield_different_carbons.png
    outpath = '../results/FBA_results/growth_rate_plots/Ppa_Yli_biomass_yield_different_carbons.png'
    if _snakemake:
        outpath = snakemake.output[0]
    fig1.savefig(outpath, dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    prepare_iLC915()
    iLC915_growth, iLC915_biomass_ylied = compute_iLC915()
    iMT1026v3_simulation_growth, iMT1026v3_simulation_biomass_yield = compute_iMT1026v3()
    iYli21_simulation_growth, iYli21_biomass_ylied = compute_yli21()
    iYali4_growth_rate, iYali4_biomass_yield = compute_yali4()
    plot_carbon_sources()

   