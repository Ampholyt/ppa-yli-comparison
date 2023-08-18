#  uses validation function from cobra and stores the results in a dataframe / csv
# takes 1.30 min
import cobra 
import os # generate storage path
import sys
import pandas as pd
import yaml

# include helper functions
sys.path.append('../scripts/')
import helperFunction as hf

config_name = 'model_config'
config_path = f'../config/{config_name}.yaml'

# load config
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# configs 
_seperator = config['seperator']
# ";" occures in error messages
_seperator = _seperator if _seperator != ';' else '\t'

def validate_models():
    ## load models
    models = hf.get_all_models() # get all models because I want to check e.coli and s.cerevisiae


    # create and store error dataframe:
    error_dict = { 'Model_Name': [],
    'SBML_FATAL': [],
    'SBML_ERROR': [],
    'SBML_SCHEMA_ERROR': [],
    'SBML_WARNING': [],
    'COBRA_FATAL': [],
    'COBRA_ERROR': [],
    'COBRA_WARNING': [],
    'COBRA_CHECK': []}

    # iterate the models and store the errors in the error_dict
    for model_name, model_path in models.items():
        # load and validate model
        model, report = cobra.io.validate_sbml_model(model_path)
        error_dict['Model_Name'].append(model_name)
        for key, value in report.items():
            error_dict[key].append(value)
    error_df = pd.DataFrame.from_dict(error_dict)

    os.makedirs('/'.join(config['results']['model_errors'].split('/')[0:-1]), exist_ok=True)  
    # save error_df
    error_df.to_csv(config['results']['model_errors'], sep=config['seperator'], index=False)

if __name__ == '__main__':
    validate_models()

