# imports + config
import cobra 
import numpy as np
import pandas as pd

import sys # append path

sys.path.append('../scripts/')
import helperFunction as hf

config = hf.load_config()

# parameters
_verbose = config['general']['verbose']
_loopless = config['general']['loopless']
_condition = 'python_script_loopless_all_glycolysis_TCA'
_snakemake = config['general']['snakemake']


_yli_maintenance = (7.8625, 7.8625)
_ppa_maintenance = (2.81, 2.81) # (according to the paper:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6838487/)



# helpful functions
def set_growth_condition(model, biomass_rxn, growth_condition, maintenance_rxn, maintenance):
    """Sets the growth condition of a given model
    @params: model (cobra model), biomass_rxn (str), growth_condition (dict), maintenance_rxn (str), maintenance (tuple of floats)"""
    # set objective function
    model.objective = biomass_rxn
    # set growth condition
    for rxn_id, flux in growth_condition.items():
        model.reactions.get_by_id(rxn_id).bounds = (-flux, flux)
    if maintenance_rxn != 'r1188':
        # set maintenance
        model.reactions.get_by_id(maintenance_rxn).bounds = maintenance

def get_reaction_fluxes(model, reactions, biomass_rxn, growth_condition, maintenance_rxn, maintenance, verbose=False):
    """Computes FBA solution of a given growth condition and returns the fluxes of the reactions of interest"""
    current_fluxes = {}

    # set growth condition
    set_growth_condition(model, biomass_rxn, growth_condition, maintenance_rxn, maintenance)
    # compute FBA solution
    # pfba solution
    pfba_solution = cobra.flux_analysis.pfba(model)
    if verbose:
        print(f'The objective value of the {model} is: {pfba_solution.fluxes[biomass_rxn]}')
    # collect pathway fluxes
    for ec_number, rxn_id in reactions.items():
        if rxn_id == "-":
            continue
        current_fluxes[ec_number] = (pfba_solution.fluxes[rxn_id],rxn_id)
    return current_fluxes

def get_fva_intervals(model, reactions, biomass_rxn, growth_condition, maintenance_rxn, maintenance, loopless=False, fva_solution="",verbose=False):
    """Computes FBA solution of a given growth condition and returns the fluxes of the reactions of interest"""
    current_fva_interval = {}
    
    # set growth condition
    set_growth_condition(model, biomass_rxn, growth_condition, maintenance_rxn, maintenance)
    
    # compute FVA solution if fva_solution == {}
    fva_solution = fva_solution if isinstance(fva_solution, pd.DataFrame) else cobra.flux_analysis.flux_variability_analysis(model, loopless=loopless)
    if verbose:
        print(f'The objective value of the {model} is: {fva_solution.loc[biomass_rxn]}')
    
    # collect pathway intervals
    for ec_number, rxn_id in reactions.items():
        if rxn_id == "-":
            continue
        current_fva_interval[ec_number] = (fva_solution.loc[rxn_id].minimum, fva_solution.loc[rxn_id].maximum)
    return fva_solution, current_fva_interval

def get_fva_flux_table(model_name, fluxes, fva_intervals, reaction_table):
    """Returns a pandas dataframe with fluxes and their corresponding fva intervals."""
    # fluxes to dataframe and key as column
    flux_df = pd.DataFrame.from_dict(fluxes, orient='index', columns=['flux', 'reaction_id'])
    # index to column with name EC
    flux_df.reset_index(inplace=True)
    flux_df.columns = ['EC', f'{model_name}_flux', f'{model_name}_reaction_id']

    # fva_intervals to dataframe
    fva_intervals = pd.DataFrame.from_dict(fva_intervals, orient='index', columns=['min', 'max'])
    # index to column with name EC
    fva_intervals.reset_index(inplace=True)
    fva_intervals.columns = ['EC', f'{model_name}_min', f'{model_name}_max']

    # left join flux_df to glycolysis table on EC
    glycolysis_result = pd.merge(pd.merge(reaction_table, flux_df, on='EC', how='left'), fva_intervals, on='EC', how='left')
    return glycolysis_result

def prepare_iLC915(iLC915_model):
    """Prepares the iLC915 model for the comparison (no carbon source)"""
    # prepare iLC915_model: needed because of internal cycles and active reactions
    # lower and upper bounds of several reactions should be set to zero according to supplementary material of iMT1026
    print('Preparing iLC915 model for comparison...')
    constraints_2_zero = ['r66','r910','r1104','r239','r111','r106','r490','r791','r243','r252','r253','r307','r308','r404','r405','r1320','r639','r640','r641','r642','r649','r650','r651','r652','r645','r646','r643','r644','r653','r654','r655','r656','r534'] 
    for rxn_id in constraints_2_zero:
        rxn = iLC915_model.reactions.get_by_id(rxn_id)
        rxn.bounds = (0,0)
    # constraint formulate uptake because of not meaningful results
    formulate_uptake = iLC915_model.reactions.get_by_id('r1143')
    formulate_uptake.bounds = (0,0)
    # important carbon sources (biotin, CO2, urea)
    important_carbon = ['r1132','r1137','r1177']

    # all carbon containing exchange reactions
    carbon_exchange_rxns = ['r1122', 'r1123', 'r1124', 'r1126', 'r1127', 'r1129', 'r1130', 'r1131', 'r1132', 'r1134', 'r1135', 'r1137', 'r1138', 'r1139', 'r1140', 'r1141', 'r1144', 'r1145', 'r1146', 'r1147', 'r1148', 'r1149', 'r1151', 'r1152', 'r1153', 'r1154', 'r1155', 'r1156', 'r1157', 'r1158', 'r1161', 'r1162', 'r1163', 'r1165', 'r1167', 'r1168', 'r1170', 'r1172', 'r1173', 'r1174', 'r1175', 'r1176', 'r1177', 'r1178']

    # set all carbon exchange reactions to 0 except important ones
    for rxn_id in carbon_exchange_rxns:
        if rxn_id not in important_carbon:
            rxn = iLC915_model.reactions.get_by_id(rxn_id)
            rxn.bounds = (0,0)

    # check definition of important rxns 
    important_rxns = ['r1133', 'r1141', 'r1145', 'r1144', 'r1148', 'r1160', 'r1137', 'r1188']

    if _verbose: 
        for rxn_id in important_rxns:
            rxn = iLC915_model.reactions.get_by_id(rxn_id)
            print(rxn.reactants, rxn.products, rxn.bounds)
            print(rxn_id, rxn.name, hf.formulaWithNames(rxn), 'rxn coefficient: ', rxn.get_coefficient(list(rxn.metabolites.keys())[0]))

    # set objective and maintenance
    biomass_rxn = 'r1339' # biomass
    # set objective: 
    iLC915_model.objective = 'r1339'

    # set maintenance reaction
    maintenance = 'r1188'
    maintenance_rxn = iLC915_model.reactions.get_by_id(maintenance)
    # maintenance_rxn.bounds = (_ppa_maintenance,_ppa_maintenance)
    # maintenance_rxn.bounds = (1,1) # according to iLC915_model
    maintenance_rxn.bounds = (1.2434,1.2434) # according to iLC915_model (1.2434 * 2.26 = 2.81)

    try: 
        iLC915_model.summary()
    except:
        print('iLC915 model: works as expected')

    gluc = iLC915_model.reactions.get_by_id('r1145')
    gluc.bounds = (-10,10)
    try: 
        iLC915_model.summary()
        print('iLC915 model: works as expected')
    except:
        print('iLC915 model: does not work as expected')
    return iLC915_model

def main():
    # model loading

    ## yli 
    ### load yali4 model
    iYali4_model = cobra.io.read_sbml_model(config['models']['yli4_corr'])
    ### load yli21 model
    iYli21_model = cobra.io.read_sbml_model(config['models']['yli21'])
    ### load iYLI647 model (no EC numbers)
    iYli647_model = cobra.io.read_sbml_model(config['models']['yli647_corr'])
    ### load yli20 model
    iYli20_model = cobra.io.read_sbml_model(config['models']['yli2.0_corr'])
    ### load iNL895 model
    iNL895_model = cobra.io.read_sbml_model(config['models']['yliNL895_corr'])
    ### load iMK735 model
    iMK735_model = cobra.io.read_sbml_model(config['models']['yliMK735_corr'])

    ## ppa
    ### load iMT1026v3 model
    iMT1026v3_model = cobra.io.read_sbml_model(config['models']['ppa1026v3'])
    ### load iLC915 model
    iLC915_model = cobra.io.read_sbml_model(config['models']['ppaiLC915'])

    # matching reactions

    ## glycolysis
    ### yli
    iYli21_glycolysis_without_branches = {
        "2.7.1.1": "R387", # id="R387" name="hexokinase (D-glucose:ATP)"
        "5.3.1.9": "R326", # id="R326" name="glucose-6-phosphate isomerase"
        "2.7.1.11": "R636", # id="R636" name="phosphofructokinase"
        "ppp->ppp": "R714", # id="R714" name="ribulose 5-phosphate 3-epimerase" # D-ribulose 5-phosphate_C5H11O8P <=> D-xylulose 5-phosphate_C5H11O8P
        "4.1.2.13": "R313", # id="R313" name="fructose-bisphosphate aldolase"
        "4.1.2.13b": "R313", # id="R256" name="D-fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase"
        "5.3.1.1": "R768", # id="R768" name="triose-phosphate isomerase"
        "1.2.1.12": "R344", # id="R344" name="glyceraldehyde-3-phosphate dehydrogenase"
        "2.7.2.3": "R642", # id="R642" name="phosphoglycerate kinase"
        "5.4.2.11": "R643", # id="R643" name="phosphoglycerate mutase"
        "4.2.1.11": "R292", # id="R292" name="enolase", 
        "2.7.1.40": "R694",# id="R694" name="pyruvate kinase"
    }

    iYali4_glycolysis_without_branches = {
        "2.7.1.1": "534", # name="hexokinase (D-glucose:ATP)" id="R_534"
        "5.3.1.9": "467", # name="glucose-6-phosphate isomerase" id="R_467"
        "2.7.1.11": "886", # name="phosphofructokinase" id="R_886"
        "4.1.2.13": "450", # name="fructose-bisphosphate aldolase" id="R_450"
        "4.1.2.13b": "450", # id="R_322" name="D-fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase" alternative to the aldolase
        "5.3.1.1": "1054", # name="triose-phosphate isomerase" id="R_1054"
        "1.2.1.12": "486", # name="glyceraldehyde-3-phosphate dehydrogenase" id="R_486"
        "2.7.2.3": "892", # name="phosphoglycerate kinase" id="R_892"
        "5.4.2.11": "893", # name="phosphoglycerate mutase" id="R_893" 
        "4.2.1.11": "366", # name="enolase" id="R_366"
        "2.7.1.40": "962",# name="pyruvate kinase" id="R_962"
    }

    iYli647_glycolysis_without_branches = {
        "2.7.1.1": "HEX1", # id="R_HEX1" name="hexokinase D glucoseATP"
        "5.3.1.9": "PGI", # R_G6PI; R_PGI id="R_PGI" name="glucose 6 phosphate isomerase"
        "2.7.1.11": "PFK", # R_PFK_3; R_PFK_2; R_PFK id="R_PFK" name="phosphofructokinase"
        "4.1.2.13": "FBA", # R_FBA id="R_FBA" name="fructose bisphosphate aldolase"
        "5.3.1.1": "TPI", # R_TPI id="R_TPI" name="triose phosphate isomerase"
        "1.2.1.12": "GAPD", # R_GAPD id="R_GAPD" name="glyceraldehyde 3 phosphate dehydrogenase"
        "2.7.2.3": "PGK", # R_PGK id="R_PGK" name="phosphoglycerate kinase"
        "5.4.2.11": "PGM", # R_PGM id="R_PGM" name="phosphoglycerate mutase"
        "4.2.1.11": "ENO", # R_ENO id="R_ENO" name="enolase"
        "2.7.1.40": "PYK",# R_PYK id="R_PYK" name="pyruvate kinase"
    }

    iMK735_glycolysis_without_branches = {
        "2.7.1.1": "HEX1", # id="R_HEX1" name="R_hexokinase__D_glucoseATP_"
        "5.3.1.9": "PGI", # id="R_PGI" name="R_glucose_6_phosphate_isomerase"
        "2.7.1.11": "PFK", # id="R_PFK" name="R_phosphofructokinase"
        "4.1.2.13": "FBA", # id="R_FBA" name="R_fructose_bisphosphate_aldolase"
        "5.3.1.1": "TPI", # id="R_TPI" name="R_triose_phosphate_isomerase"
        "1.2.1.12": "GAPD", # id="R_GAPD" name="R_glyceraldehyde_3_phosphate_dehydrogenase"
        "2.7.2.3": "PGK", # id="R_PGK" name="R_phosphoglycerate_kinase"
        "5.4.2.11": "PGM", # id="R_PGM" name="R_phosphoglycerate_mutase" 
        "4.2.1.11": "ENO", # id="R_ENO" name="R_enolase"
        "2.7.1.40": "PYK",# id="R_PYK" name="R_pyruvate_kinase"
    }

    iNL865_glycolysis_without_branches = {
        "2.7.1.1": "r_0573", # id="r_0573" name="hexokinase (D-glucose:ATP)"
        "5.3.1.9": "r_0505", # id="r_0505" name="glucose-6-phosphate isomerase"
        "2.7.1.11": "r_0859", # id="r_0859" name="phosphofructokinase" sboTerm="SBO:0000176"
        "4.1.2.13": "r_0484", # id="r_0484" name="fructose-bisphosphate aldolase" sboTerm="SBO:0000176"
        "5.3.1.1": "r_1041", # id="r_1041" name="triose-phosphate isomerase" sboTerm="SBO:0000176"
        "1.2.1.12": "r_0525", # id="r_0525" name="glyceraldehyde-3-phosphate dehydrogenase" sboTerm="SBO:0000176"
        "2.7.2.3": "r_0865", # id="r_0865" name="phosphoglycerate kinase" sboTerm="SBO:0000176"
        "5.4.2.11": "r_0866", # id="r_0866" name="phosphoglycerate mutase" sboTerm="SBO:0000176"
        "4.2.1.11": "r_0398", # id="r_0398" name="enolase" sboTerm="SBO:0000176"
        "2.7.1.40": "r_0941",# id="r_0941" name="pyruvate kinase" sboTerm="SBO:0000176"
    }

    iYL_20_glycolysis_without_branches = {
        "2.7.1.1": "R0371", # id="R0371" name="ATP:beta-D-glucose 6-phosphotransferase"
        "5.3.1.9": "R0376", # id="R0376" name="alpha-D-Glucose 6-phosphate ketol-isomerase"
        "2.7.1.11": "R0379", # id="R0379" name="ATP:D-fructose-6-phosphate 1-phosphotransferase"
        "4.1.2.13": "R0368", # id="R0368" name="beta-D-fructose-1,6-bisphosphate D-glyceraldehyde-3-phosphate-lyase (glycerone-phosphate-forming)"
        "5.3.1.1": "R0366", # id="R0366" name="D-glyceraldehyde-3-phosphate aldose-ketose-isomerase"
        "1.2.1.12": "R0367", # id="R0367" name="D-glyceraldehyde-3-phosphate:NAD+ oxidoreductase (phosphorylating)"
        "2.7.2.3": "R0369", # id="R0369" name="ATP:3-phospho-D-glycerate 1-phosphotransferase"
        "5.4.2.11": "R0370", # id="R0370" name="2-Phospho-D-glycerate 2,3-phosphomutase"
        "4.2.1.11": "R0360", # id="R0360" name="2-phospho-D-glycerate hydro-lyase (phosphoenolpyruvate-forming)"
        "2.7.1.40": "R0382",# id="R0382" name="ATP:pyruvate 2-O-phosphotransferase"
    }

    ### ppa 

    iMT1026v3_glycolysis_without_branches = {
        "2.7.1.1": "HEX1", # id="R_HEX1" name="hexokinase (D-glucose:ATP)"
        "5.3.1.9": "PGI", # id="R_PGI" name="glucose-6-phosphate isomerase"
        "5.1.3.3(a-D-g->b-D-g6p)": "GLUK", # +2.7.1.1 a-D-g -> b-D-g6p id="R_GLUK" Alternative pathway from glucose to fructose-6-phosphate
        "5.1.3.9(b-D-g6P<=>b-D-f6P)": "G6PI3", # id="R_G6PI3" name="glucose-6-phosphate isomerase (b-D-glucose-6P <=> b-D-fructose-6P)"
        "5.1.3.9(g6p_c<=>g6p_B_c)": "G6PI", # id="R_G6PI" name="glucose-6-phosphate isomerase (g6p_c <=> g6p_B_c)"
        "1.1.1.49(gly->ppp)": "G6PDH2", # id="G6PDH2" name="glucose 6-phosphate dehydrogenase" g6p -> PPP
        "5.1.3.1(ru5p_D_c<=>xu5p_D_c)": "RPE", # id="RPE" name="ribulose 5-phosphate 3-epimerase" ru5p_D_c <=> xu5p_D_c
        "2.2.1.1(ppp->gly)": "TKT2", # id="TKT2" name="transketolase" PPP from ppp -> f6p + g3p
        "2.7.1.11": "PFK", # id="R_PFK" name="phosphofructokinase"
        "4.1.2.13": "FBA", # id="R_FBA" name="fructose-bisphosphate aldolase"
        "4.1.2.13b": "FBA2", # id="R_FBA2" name="D-Fructose 1-phosphate D-glyceraldehyde-3-phosphate-lyase" (alternative zu FBA)
        "5.3.1.1": "TPI", # id="R_TPI" name="triose-phosphate isomerase"
        "1.2.1.12": "GAPD", # id="R_GAPD" name="glyceraldehyde-3-phosphate dehydrogenase"
        "2.7.2.3": "PGK", # id="R_PGK" name="phosphoglycerate kinase"
        "5.4.2.11": "PGM", # id="R_PGM" name="phosphoglycerate mutase" # iMT1026v3: EC: 5.4.2.1
        "4.2.1.11": "ENO", # id="R_ENO" name="enolase"
        "2.7.1.40": "PYK",# id="R_PYK" name="pyruvate kinase" 
    }

    iLC915_glycolysis_without_branches = {
        "2.7.1.1": "r552", # id="R_r552" name="ATP:alpha-D-glucose 6-phosphotransferase"
        "5.3.1.9": "r984", # id="R_r984" name="beta-D-Glucose 6-phosphate ketol-isomerase"
        "PPP->glycer3p": "r319", # id="R_r319" name="sedoheptulose-7-phosphate:D-glyceraldehyde-3-phosphate"
        "2.7.1.11": "r557", # id="R_r557" name="ATP:D-fructose-6-phosphate 1-phosphotransferase"
        "4.1.2.13": "r872", # id="R_r872" name="beta-D-fructose-1,6-bisphosphate D-glyceraldehyde-3-phosphate-lyase" 
        "5.3.1.1": "r978", # id="R_r978" name="D-glyceraldehyde-3-phosphate aldose-ketose-isomerase"
        "1.2.1.12": "r146", # id="R_r146" name="D-glyceraldehyde-3-phosphate:NAD+ oxidoreductase (phosphorylating)" 
        "2.7.2.3": "r596", # id="R_r596" name="ATP:3-phospho-D-glycerate 1-phosphotransferase"
        "5.4.2.11": "r988", # id="R_r988" name="2-Phospho-D-glycerate 2,3-phosphomutase" 
        "4.2.1.11": "r891", # id="R_r891" name="2-phospho-D-glycerate hydro-lyase (phosphoenolpyruvate-forming)" 
        "2.7.1.40": "r576",# id="R_r576" name="ATP:pyruvate 2-O-phosphotransferase"
    }

    ## citrate cycle

    ### ppa
    # iMT1026v3
    iMT1026v3_citrate_cycle = {
    '6.4.1.1': 'PC', # NOT MITOCHORIAL pyruvate carboxylase (Pyruvate -> Oxaloacetate)
    '1.1.1.37': 'MDHm', # malate dehydrogenase, mitochondrial (Malate -> Oxaloacetate)
    'pyruvate_mitochodrial_transport': 'PYRt2m', # pyruvate mitochondrial transport via proton symport
    '1.2.4.1': 'PDHa1', # pyruvate dehydrogenase (lipoamide)
    '1.2.4.1_2': 'PDHa2', # pyruvate dehydrogenase (lipoamide)
    '2.3.1.12': 'PDHbrm', # pyruvate dehydrogenase (dihydrolipoamide) reversible (mitochondrial) 
    '2.3.3.1': 'CSm', # citrate synthase (mitochondrial)
    '4.2.1.3': 'ACONTm', # aconitase (mitochondrial)
    '1.1.1.41': 'ICDHxm', # isocitrate dehydrogenase (NAD+)
    '1.2.4.2_1': 'AKGDH1', # 2-oxoglutarate dehydrogenase
    '1.2.4.2_2': 'AKGDH2', # 2-oxoglutarate dehydrogenase
    '2.3.1.61': 'AKGDbm', # oxoglutarate dehydrogenase (dihydrolipoamide S-succinyltransferase) (mitochondrial)
    '6.2.1.4': 'SUCOASm', # Succinate--CoA ligase (ADP-forming) # according to KEGG: 6.2.1.4: GDP-forming
    '6.2.1.5': 'ITCOALm', # Itaconate--CoA ligase (ADP-forming), mitochondrial id="R_ITCOALm" (?)
    '1.3.5.1': 'SUCD2_u6m', # succinate dehydrogenase (ubiquinone-6), mitochondrial
    '4.2.1.2': 'FUMm', # fumarase (mitochondrial)
    }

    ### yli
    # iYli21
    iYli21_citrate_cycle = {
    '6.4.1.1': 'R690', # NOT MITOCHONRIAL pyruvate carboxylase (Pyruvate -> Oxaloacetate)
    'pyruvate_mitochodrial_transport': 'R1311', # pyruvate mitochondrial transport via proton symport
    '1.1.1.37': 'R533', # (m) malate dehydrogenase (NAD+) (Malate -> Oxaloacetate)
    '1.2.4.1': 'R693', # (m) pyruvate dehydrogenase (Pyruvate -> Acetyl-CoA)
    '2.3.3.1': 'R236', # (m) citrate synthase (Oxaloacetate + Acetyl-CoA -> Citrate)
    '4.2.1.3': 'R238', # (m) aconitase (citrate to cis-aconitate(3-)) # according to Expasy EC = 4.2.1.3
    '1.1.1.41': 'R487', # (m) isocitrate dehydrogenase (NAD+)
    '1.1.1.42': 'R1385', # (m) isocitrate dehydrogenase (NAD+)
    '1.2.4.2': 'R614', # (m) oxoglutarate dehydrogenase (lipoamide)
    '2.3.1.61': 'R613', # (m) oxoglutarate dehydrogenase (dihydrolipoamide S-succinyltransferase)
    '6.2.1.4': 'R741', # (m) succinate-CoA ligase (ADP-forming)
    '6.2.1.5': 'R496', # (m) itaconate-CoA ligase (ADP-forming) id="R496" (?)
    '1.3.5.1': 'R740', # (m) succinate dehydrogenase (ubiquinone-6)
    '4.2.1.2': 'R314', # (m) fumarase (fumarate to malate)
    '1.1.1.37': 'R533', # (m) malate dehydrogenase (malate to oxaloacetate)
    }

    # iYali4
    iYali4_citrate_cycle = {
    '6.4.1.1': '958', # pyruvate carboxylase (id:R_958) (Pyruvate -> Oxaloacetate)
    '1.1.1.37': '713', # malate dehydrogenase (id:R_713) (Malate -> Oxaloacetate)
    'pyruvate_mitochodrial_transport': '2034', # name="pyruvate transport" id="R_2034"
    '1.2.4.1': '961', # name="pyruvate dehydrogenase" id="R_961"
    '1.2.4.1_2': '-', # pyruvate dehydrogenase (lipoamide)
    '2.3.1.12': '-', # pyruvate dehydrogenase (dihydrolipoamide) reversible (mitochondrial) 
    '2.3.3.1': '300', # citrate synthase (mitochondrial) id="R_300"
    '4.2.1.3': '2305', # aconitase (mitochondrial) name="cis-aconitate(3-) to isocitrate" id="R_2305" (?)
    '1.1.1.41': '658', # isocitrate dehydrogenase (NAD+) id="R_658"
    '1.2.4.2': '832', # name="oxoglutarate dehydrogenase (lipoamide)" id="R_832"
    '1.2.4.2_2': '-', # 
    '2.3.1.61': '831', # oxoglutarate dehydrogenase (dihydrolipoamide S-succinyltransferase) id="R_831"
    '6.2.1.4': '1022', # succinate-CoA ligase (ADP-forming)" id="R_1022" 
    '6.2.1.5': '668', # name="itaconate-CoA ligase (ADP-forming)" id="R_668" (?)
    '1.3.5.1': '1021', # succinate dehydrogenase (ubiquinone-6) (m) id="R_1021"
    '4.2.1.2': '451', # fumarase (mitochondrial) id="R_451"
    }

    # iYLI647
    iYli647_citrate_cycle = {
    '6.4.1.1': 'PC', # (id:R_PC) pyruvate carboxylase (Pyruvate -> Oxaloacetate)
    '1.1.1.37': 'MDHm', # (id:R_MDHm) malate dehydrogenase, mitochondrial (Malate -> Oxaloacetate)  
    'pyruvate_mitochodrial_transport': 'PYRt2m', # id="R_PYRt2m" name="pyruvate mitochondrial transport via proton symport" Pyruvate c -> pyruvate mito
    '1.2.4.1': 'PDHm', # pyruvate dehydrogenase id="R_PDHm"
    '1.2.4.1_2': '-', # pyruvate dehydrogenase (lipoamide)
    '2.3.1.12': '-', # pyruvate dehydrogenase (dihydrolipoamide) reversible (mitochondrial) 
    '2.3.3.1': 'CSm', # citrate synthase (mitochondrial) id="R_CSm"
    '4.2.1.3': 'ACONT', # aconitase (mitochondrial) id="R_ACONT"
    '1.1.1.41': 'ICDHxm', # isocitrate dehydrogenase (NAD+) id="R_ICDHxm"
    '1.2.4.2_1': 'AKGDam', # 2-oxoglutarate dehydrogenase id="R_AKGDam" name="oxoglutarate dehydrogenase lipoamide"
    '1.2.4.2_2': '-', # 2-oxoglutarate dehydrogenase 
    '2.3.1.61': 'AKGDbm', # oxoglutarate dehydrogenase dihydrolipoamide S succinyltransferase id="R_AKGDbm" (mitochondrial)
    '6.2.1.4': 'SUCOASm', # Succinate--CoA ligase (GDP-forming)
    '6.2.1.5': 'ITCOALm', # Itaconate CoA ligase ADP forming mitochondrial id="R_ITCOALm" (?)
    '1.3.5.1': 'SUCD2_u6m', # succinate dehydrogenase ubiquinone 6 mitochondrial id="R_SUCD3_u6m"
    '4.2.1.2': 'FUMm', # fumarase mitochondrial id="R_FUMm"
    }

    # combining glycolysis and citrate cycle
    ## ppa (only iMT1026v3 and glycolysis of iLC915)
    iMT1026v3_combined = {**iMT1026v3_glycolysis_without_branches, **iMT1026v3_citrate_cycle}

    ## yli
    yali4_combined = {**iYali4_glycolysis_without_branches, **iYali4_citrate_cycle}
    iYli21_combined = {**iYli21_glycolysis_without_branches, **iYli21_citrate_cycle}
    yli647_combined = {**iYli647_glycolysis_without_branches, **iYli647_citrate_cycle}

    # reading reaction table
    reaction_table = pd.read_csv(config['experiments']['reaction_table'], sep=config['seperator'])

    # setting up the models and getting the fva results

    ## ppa

    # # iLC915 results
    # model_name = 'iLC915'
    # print(f"{model_name}: \n")
    # biomass_rxn = 'r1339' # r1339 Growth # r1187 biomass formation
    # maintenance_rxn = 'r1188' # Maintenance
    # growth_condition = {
    #     'r1145': 100, # glucose
    #     'r1144': 0, # fructose
    #     'r1148': 0, # glycerol
    #     'r1160': 1000, # O2
    # }
    # iLC915_model = prepare_iLC915(iLC915_model)


    # # fluxes and fva intervals
    # current_fluxes = get_reaction_fluxes(iLC915_model, iLC915_glycolysis_without_branches, biomass_rxn, growth_condition, maintenance_rxn, _ppa_maintenance)
    # iLC915_fva_solution, current_fva_intervals = get_fva_intervals(iLC915_model, iLC915_glycolysis_without_branches, biomass_rxn, growth_condition, maintenance_rxn, _ppa_maintenance, _loopless, "", _verbose)
    # # get flux table
    # reaction_table = get_fva_flux_table(model_name, current_fluxes, current_fva_intervals, reaction_table)
    # # output
    # # reaction_table.to_csv(f"{config['results']['metabolic_fluxes']}/{_condition}_flux_fva.csv", sep=config['seperator'], index=False)

    # # iMT1026v3 results
    model_name = 'iMT1026v3'
    print(f"{model_name}: \n")
    biomass_rxn = 'growth' # BIOMASS_glyc # (id:BIOMASS) Biomass composition # growth biomass_c --> biomass_e
    maintenance_rxn = 'ATPM' # 0.55 
    growth_condition = {
        'Ex_glc_D': 100, # glucose
        'Ex_fru': 0, # fructose
        'Ex_glyc': 0, # glycerol
        'Ex_o2': 1000, # O2
    }

    # fluxes and fva intervals
    current_fluxes = get_reaction_fluxes(iMT1026v3_model, iMT1026v3_combined, biomass_rxn, growth_condition, maintenance_rxn, _ppa_maintenance)
    iMT1026v3_fva_solution, current_fva_intervals = get_fva_intervals(iMT1026v3_model, iMT1026v3_combined, biomass_rxn, growth_condition, maintenance_rxn, _ppa_maintenance, _loopless, "", _verbose)
    # get flux table
    reaction_table = get_fva_flux_table(model_name, current_fluxes, current_fva_intervals, reaction_table)
    # output
    # reaction_table.to_csv(f"{config['results']['metabolic_fluxes']}/{_condition}_flux_fva.csv", sep=config['seperator'], index=False)
    # reaction_table.to_csv(config['results']['metabolic_comparison_condition'], sep=config['seperator'], index=False)


    ## yli

    ### iYli647 results
    model_name = "iYLI647"
    print(f"{model_name}: \n")
    biomass_rxn = 'biomass_C'
    maintenance_rxn = 'ATPM'
    growth_condition = {
        'EX_glc(e)': 100, # glucose
        'EX_fru(e)': 0, # fructose
        'EX_glyc(e)': 0, # glycerol
    }

    # For glycolysis and citrate cycle
    iYli647_glycolysis_and_citrate_cycle = {**iYli647_glycolysis_without_branches, **iYli647_citrate_cycle}

    current_fluxes = get_reaction_fluxes(iYli647_model, iYli647_glycolysis_and_citrate_cycle, biomass_rxn, growth_condition, maintenance_rxn, _yli_maintenance)
    iYli647_fva_solution, current_fva_intervals = get_fva_intervals(iYli647_model, iYli647_glycolysis_and_citrate_cycle, biomass_rxn, growth_condition, maintenance_rxn, _yli_maintenance, _loopless, "", _verbose)
    # add both to the csv table
    reaction_table = get_fva_flux_table(model_name, current_fluxes, current_fva_intervals, reaction_table)
    # reaction_table.to_csv(f"{config['results']['metabolic_fluxes']}/{_condition}_flux_fva.csv", sep=config['seperator'], index=False)
    # reaction_table.to_csv(config['results']['metabolic_comparison_condition'], sep=config['seperator'], index=False)


    ### iYali4 results
    model_name = 'iYali4'
    print(f"{model_name}: \n")
    biomass_rxn = 'biomass_C'
    maintenance_rxn = 'xMAINTENANCE'
    growth_condition = {
        '1714': 100, # glucose
        '1709': 0, # fructose
        '1808': 0, # glycerol
    }

    current_fluxes = get_reaction_fluxes(iYali4_model, yali4_combined, biomass_rxn, growth_condition, maintenance_rxn, _yli_maintenance)
    iYali4_fva_solution, current_fva_intervals = get_fva_intervals(iYali4_model, yali4_combined, biomass_rxn, growth_condition, maintenance_rxn, _yli_maintenance, _loopless, "", _verbose)
    # add both to the csv table
    reaction_table = get_fva_flux_table(model_name, current_fluxes, current_fva_intervals, reaction_table)
    # reaction_table.to_csv(f"{config['results']['metabolic_fluxes']}/{_condition}_flux_fva.csv", sep=config['seperator'], index=False)
    # reaction_table.to_csv(config['results']['metabolic_comparison_condition'], sep=config['seperator'], index=False)


    ### iYli21 results
    model_name = 'iYli21'
    print(f"{model_name}: \n")
    biomass_rxn = 'biomass_C'
    maintenance_rxn = 'xMAINTENANCE'
    growth_condition = {
        'R1070': 100, # glucose
        'R1065': 0, # fructose
        'R1141': 0, # glycerol
    }
    current_fluxes = get_reaction_fluxes(iYli21_model, iYli21_combined, biomass_rxn, growth_condition, maintenance_rxn, _yli_maintenance)
    iYli21_fva_solution, current_fva_intervals = get_fva_intervals(iYli21_model, iYli21_combined, biomass_rxn, growth_condition, maintenance_rxn, _yli_maintenance, _loopless, "", _verbose)
    # add both to the csv table
    reaction_table = get_fva_flux_table(model_name, current_fluxes, current_fva_intervals, reaction_table)
    # reaction_table.to_csv(f"{config['results']['metabolic_fluxes']}/{_condition}_flux_fva.csv", sep=config['seperator'], index=False)
    reaction_table.to_csv(config['results']['metabolic_comparison_condition'], sep=config['seperator'], index=False)

    print('finished metabolic analysis and stored the Flux Variability Analysis (FVA) results in a csv file.')

if __name__ == "__main__":
    main()


