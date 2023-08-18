## growth plot yli (FBA)
rule yli_exp_data:
    input:
        config['results']['model_errors'],
        config['results']['model_numbers'],
    output:
        'results/FBA_results/yli_glucose_growth_comparison.png',
        'results/FBA_results/{condition_name}_simulated_growth_rates.csv'.format(condition_name = config['general']['condition_name']),
    log:
        "results/logs/yli_exp_data.log"
    conda:
        "../envs/fba_results.yaml"
    script:
        '../scripts/yli_fba_plot.py'

## RSS plot of yli
rule RSS_plot_yli:
    output: 
        'results/FBA_results/residual_sum_of_squares/top3_RSS_model_comparison.png'
    log:
        "results/logs/RSS_plot_yli.log"
    conda:
        "../envs/quality_control.yaml"
    script:
        "../scripts/top3_rss.py"

## growth plot ppa (FBA)
rule ppa_exp_data_281:
    output:
        'results/FBA_results/ppa_glucose_growth_comparison_281.png',
        'results/FBA_results/ppa_CO2_production_281.png',
    params:
        2.81,
    log:
        "results/logs/ppa_exp_data_281.log"
    conda:
        "../envs/fba_results.yaml"
    script: 
        '../scripts/ppa_fba_plot.py'

rule ppa_exp_data_226:
    output:
        'results/FBA_results/ppa_glucose_growth_comparison_226.png',
        'results/FBA_results/ppa_CO2_production_226.png',
    params:
        2.26,
    log:
        "results/logs/ppa_exp_data_226.log"
    conda:
        "../envs/fba_results.yaml"
    script: 
        '../scripts/ppa_fba_plot.py'

## growth comparison on different carbon sources (FBA)
rule carbon_comparison:
    output:
        'results/FBA_results/carbon_comparison.png',
    log:
        "results/logs/carbon_comparison.log"
    conda:
        "../envs/fba_results.yaml"
    script:
        '../scripts/carbon_comparison.py'
