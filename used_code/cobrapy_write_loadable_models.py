# This script reads the loadable models and writes them to a seperate folder. This should fix some of the occurring warnings.
# Gets the loadable models from the helperFunction.py script.
import cobra 
import os
import sys 

sys.path.append('../scripts/')
import helperFunction as hf

config = hf.load_config()

def read_and_write_model(input_path, outdir=None):
    split_path = input_path.split('/')
    file_name = split_path[-1].split('.')[0]
    model = cobra.io.read_sbml_model(input_path)
    outpath = os.path.join('/'.join(split_path[:-1]), file_name + '_read_write.xml') if not outdir else outdir + file_name + '_read_write.xml'
    os.makedirs(os.path.dirname(outpath), exist_ok=True) # Create the directory if it does not exist
    cobra.io.write_sbml_model(model,outpath)
    return outpath

def main():
    loadable_models = hf.get_loadable_models()

    for model_name, model_path in loadable_models.items():
        print(model_name)
        read_and_write_model(model_path, config['models']['cobra_wrote']['directory'])
        print('-'*50)

if __name__ == '__main__':
    main()