# Quality control of the models
rule install_memote:
    output:
        "results/memote/installed.txt"
    log:
        "results/logs/memote_installation.log"
    params:
        output_dir = "results/memote",
    shell:
        "pip install memote && mkdir -p {params.output_dir} && echo memote was installed successfully > {output}"

rule produce_memote_report:
    input:
        "results/memote/installed.txt"
    output:
        "results/memote/{model_name}_report.html"
    params:
        model_name=lambda wildcards: wildcards.model_name,
        model_path=lambda wildcards: all_models[wildcards.model_name],
    log:
        "results/logs/memote_report_{model_name}.log"
    shell:
        "memote report snapshot --filename {output} {params.model_path}"
    ##### Error with iMT1026Chan2017_report.html and yli647_report.html
    ### iMT1026Chan2017_report.html
    # error: Reference: L3V1 Fbc V3 Section 3.4
    # error:  Encountered '_' when expecting a capital letter. The chemicalFormula 'Apo_ficytcCharge2' has incorrect syntax.
    # error:  - Category: SBML component consistency, Severity: 2

    ### yli647_report.html
    # ../../../../../../miniconda3/envs/sn/lib/python3.9/site-packages/memote/suite/tests/test_consistency.py FFFatal Python error: Aborted
    # Current thread 0x00000002031322c0 (most recent call first):
    # File "/Users/ampholyt/miniconda3/envs/sn/lib/python3.9/site-packages/swiglpk/swiglpk.py", line 569 in glp_intopt

rule run_model_validation:
    output:
        config['results']['model_errors']
    log:
        "results/logs/model_validation.log"
    conda:
        "../envs/quality_control.yaml"
    script:
        "../scripts/run_model_validation.py"

rule generate_model_number_table:
    output:
        config['results']['model_numbers']
    log:
        "results/logs/generate_model_number.log"
    conda:
        "../envs/quality_control.yaml"
    script:
        "../scripts/generate_model_numbers_table.py"

