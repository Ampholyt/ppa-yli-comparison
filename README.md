# Comparing genome-scale models of the non-conventional yeast <em>Y. lipolytica</em> and <em>P. pastoris</em>	
## 1. Installations
### Installing default environment
1. installing conda 
2. installing mamba
```bash
conda create -n carbon -c conda-forge python=3.8 mamba
conda activate carbon
mamba env update -n carbon --file workflow/envs/default_env.yaml
```

### Installing memote (if this task fails in the workflow because of your local environment)
1. If needed create new envinroment (e.g. virtual python environment or conda/mamba)


1.1 generate virtual python environment with cobrapy as name
```bash
python3 -m venv cobrapy
```
2. Install memote using pip
2.1 
```bash 
pip install memote
```
3. Use memote to analyse sbml model
```bash 
memote report snapshot <path to xml-model>
```

### Supplementary data
Here the presented and cleaned version of the produced data is stored. One can find all the images and tables used to produce the thesis. The other directories contain either code or the same data in less clean form e.g. the metabolic map was created manually and is therefore only available in the ~/supplementary_data directory.

### Starting the analysis (generating the results for the analysis)
0. Modify the config file (~/workflow/config/model_config.yaml)
1. Start snakemake with n cores (in the directory of the Snakefile)
```bash
snakemake --cores <nuber_of_threads> --keep-going --use-conda
```
- keep-going: this flag is required because some of the provided models will lead to errors which will stop the workflow if this flag is not set
- use-conda: allows the workflow to use conda as a package management tool

### Investigating the comparison of the central carbon metabolism 
The notebook provided in the ~/workflow/scripts directory can be run and depending on your preference with or without loopless flag (see config) it might take over 30 minutes (depending on your hardware).

### Investigate the analysis 
1. Check the results directory
2. For the code required for the analysis see ~/workflow/scripts and for additional and not required code see ~/used_code directory.



