# script needed for wrong translated models (translation between .mat and .xml type integrated errors if there exist multiple underscores in the model name)
# And the error was not visible in jupyter (kernel died)
import cobra
import yaml


# config
config_name = 'model_config'
config_path = f'../config/{config_name}.yaml'

# load config
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)


# load model
iyli647_model_path = config['models']['yli647_test']
iyli647_model = cobra.io.read_sbml_model(iyli647_model_path)

# default growth:
sol = iyli647_model.optimize()
print(sol.objective_value)