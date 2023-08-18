### EFM (upset plot, number comparison)
rule efm_decomp_iLC915:
    output:
        'results/EFM_decomp_iLC915/scaled_coefs.png',
        'results/EFM_decomp_iLC915/active_reactions.png',
        'results/EFM_decomp_iLC915/upset_plot_EFM1-4.png',
    log:
        'results/logs/EFM_decomp_iLC915.log'
    conda:
        "../envs/efm_decomp.yaml"
    script:
        "../scripts/efm_decomp_iLC915.py"
