rule fastqc:
    input:
        "data/reads/{library}.fq.gz"
    output:
        html="results/LCMSeq/quality_control/fastqc_report/{library}_fastqc.html",
        zip="results/LCMSeq/quality_control/fastqc_report/{library}_fastqc.zip"
    params: "--quiet"
    log:
        config["Log_Dir"] + "/fastqc/{library}.log"
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    script:
        "../wrappers/fastqc.py"

rule multiqc:
    input:
        expand("results/LCMSeq/quality_control/fastqc_report/{sample}_{end}_fastqc.html",
               sample=config["Samples"], end=["R1", "R2"])
    output:
        "results/LCMSeq/quality_control/fastqc_report/multiqc.html"
    log:
        config["Log_Dir"] + "/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    script:
        "../wrappers/multiqc.py"