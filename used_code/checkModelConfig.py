# This skrip should test the model configuration: iterate over the models and verify that the config's path leads to a file.
import os
import yaml

config_name = 'model_config'

def checkModels(config):
    """Iterates the models in the config file and checks if the path exists."""
    model_path = config['models']
    for model, path in model_path.items():
        print(f'Checking model {model}...')
        if os.path.exists(path):
            print(f'Path exists.')
        else:
            print(f'Path does not exist. \nThe configuration file is not set correctly. Please check the path of {model} there is no model at: {path}.\nThe configuration file is at: {config_path}')
            return False
    return True

# main
if __name__ == '__main__':
    # config
    config_path = f'../config/{config_name}.yaml'

    # load config
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    # check models
    if (checkModels(config)):
        print('All models exist.')
    else:
        raise Exception('Some models do not exist. Check the output above.')