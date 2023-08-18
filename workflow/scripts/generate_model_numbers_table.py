import cobra
import sys  
import os # generate path  
import pandas as pd
import yaml

# include helper functions
sys.path.append('../scripts/')
import helperFunction as hf

# get config
config = hf.load_config()

def create_model_numbers_table(models, outpath, seperator):
    '''Uses models dictionary to create one dataframe with the model numbers'''
    numberDFs = []
    # iterate all models 
    for model_name, model_path in models.items():
        # load and validate model
        model, report = cobra.io.validate_sbml_model(model_path)
        # check if model is nontype 
        if isinstance(model,type(None)):
            print(f'{model_name} is NoneType')
            continue
        numberDFs.append(hf.numberDf(model,model_name, config['results']['model_numbers_directory'], seperator, True, False))

    # concat all dataframes
    modelNumbers = pd.concat(numberDFs)

    # save modelNumbers
    modelNumbers.to_csv(outpath, sep=seperator, index=False)
    return modelNumbers

if __name__ == '__main__':
    seperator = config['seperator'] # ';' is fine

    ## load models
    models = hf.get_all_models()
    
    # generate model numbers table
    modelNumbers = create_model_numbers_table(models, config['results']['model_numbers'], seperator)