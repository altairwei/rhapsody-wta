## RSEQC


rule rseqc_gtf2bed:
    input:
        config["Reference_Genome"]["Annotation"],
    output:
        bed="results/LCMSeq/align/rseqc/annotation.bed",
        db=temp("results/LCMSeq/align/rseqc/annotation.db"),
    log:
        config["Log_Dir"] + "/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bed="results/LCMSeq/align/rseqc/annotation.bed",
    output:
        "results/LCMSeq/align/rseqc/{sample}.junctionanno.junction.bed",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_junction_annotation/{sample}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bed="results/LCMSeq/align/rseqc/annotation.bed",
    output:
        "results/LCMSeq/align/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/LCMSeq/align/rseqc/{sample}.stats.txt",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_stat/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bed="results/LCMSeq/align/rseqc/annotation.bed",
    output:
        "results/LCMSeq/align/rseqc/{sample}.infer_experiment.txt",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_infer/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bed="results/LCMSeq/align/rseqc/annotation.bed",
    output:
        "results/LCMSeq/align/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_innerdis/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bed="results/LCMSeq/align/rseqc/annotation.bed",
    output:
        "results/LCMSeq/align/rseqc/{sample}.readdistribution.txt",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_readdis/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/LCMSeq/align/rseqc/{sample}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_readdup/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/LCMSeq/align/rseqc/{sample}.readgc.GC_plot.pdf",
    priority: 1
    log:
        config["Log_Dir"] + "/rseqc/rseqc_readgc/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule align_multiqc:
    input:
        expand(
            "results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.junctionanno.junction.bed",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.infer_experiment.txt",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.stats.txt",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.readdistribution.txt",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.readdup.DupRate_plot.pdf",
            sample=config["Samples"],
        ),
        expand(
            "results/LCMSeq/align/rseqc/{sample}.readgc.GC_plot.pdf",
            sample=config["Samples"],
        ),
        expand(
            config["Log_Dir"] + "/rseqc/rseqc_junction_annotation/{sample}.log",
            sample=config["Samples"],
        ),
    output:
        "results/LCMSeq/align/rseqc/multiqc_report.html",
    log:
        config["Log_Dir"] + "/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    script:
        "../wrappers/multiqc.py"
