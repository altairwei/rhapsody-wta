rule trimmomatic_pe:
    input:
        r1="data/reads/{sample}_R1.fq.gz",
        r2="data/reads/{sample}_R2.fq.gz"
    output:
        r1="results/LCMSeq/trimmed_reads/{sample}_R1.fq.gz",
        r2="results/LCMSeq/trimmed_reads/{sample}_R2.fq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="results/LCMSeq/trimmed_reads/{sample}_R1_unpaired.fq.gz",
        r2_unpaired="results/LCMSeq/trimmed_reads/{sample}_R2_unpaired.fq.gz"
    log: "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=[
            "ILLUMINACLIP:data/misc/TruSeq3-PE-2.fa:2:30:10:2:true",
            "SLIDINGWINDOW:4:20",
            "TRAILING:20",
            "MINLEN:20"
        ],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        32
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    conda:
        "../envs/trimmomatic.yaml"
    script:
        "../wrappers/trimmomatic_pe.py"

rule fastqc_trimmed:
    input:
        "results/LCMSeq/trimmed_reads/{library}.fq.gz"
    output:
        html="results/LCMSeq/quality_control/trimmed_fastqc_report/{library}_fastqc.html",
        zip="results/LCMSeq/quality_control/trimmed_fastqc_report/{library}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{library}.log"
    threads: 1
    conda:
        "../envs/fastqc.yaml"
    script:
        "../wrappers/fastqc.py"

rule multiqc_trimmed:
    input:
        expand("results/LCMSeq/quality_control/trimmed_fastqc_report/{sample}_{end}_fastqc.html",
               sample=config["Samples"], end=["R1", "R2"])
    output:
        "results/LCMSeq/quality_control/trimmed_fastqc_report/multiqc.html"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    script:
        "../wrappers/multiqc.py"