## Snakemake workflow for the results of the analysis and comparison of the two non-conventional yeasts P. pastoris and Y. lipolytica
This workflow will generate all images used in the thesis "Comparing genome-scale models of the non-conventional yeast Y. lipolytica and P. pastoris"

### Getting started
0. Install the package manager [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
1. Get your environment ready and install the default conda environment with the name conp_ppa_yli (env/default_env.yaml)
```bash
conda create -n carbon -c conda-forge python=3.8 mamba
conda activate carbon
mamba env update -n carbon --file workflow/envs/default_env.yaml
```
2. Run the workflow
```bash
snakemake --cores <nuber_of_threads> --keep-going --use-conda
```
- keep-going: this flag is required because some of the provided models will lead to errors which will stop the workflow if this flag is not set
- use-conda: allows the workflow to use conda as package management tool

This workflow is testes with Python 3.8

### Workflow: 
- Validate the models -> error table (run_model_validation.py)
- Investigate predefined growth condition for all the loadable models -> growth result table
- Species specific model analysis
    - Yli:
        - Experimental results -> Image (yli_fba_plot.py)
        - RSS -> best performing models 
        - (Currency analysis)
    - Ppa:
        - Experimental results (Ppa) -> Image (ppa_fba_plot.py)
        - (RSS) (not necessary because both models will be used for further analysis)
        - (Currency analysis) 
    
- Species Comparison:
    - Growth on different growth conditions -> Image (biomass yield production per carbon source in C mol/gDCW/h) (carbon_comparison.py)
    - Flux Variability -> Table of metabolic flux ranges (turned into map manually)

    - EFM decomposition analysis (iLC915)
        - analysis of numbers of active reactions
        - table of grouped EFMs: proportion to biomass, oxygen uptake and carbon uptake
