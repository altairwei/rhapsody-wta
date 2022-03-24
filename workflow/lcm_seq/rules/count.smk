###############################################################################
# TRANSCRIPTOME QUATIFICATION
###############################################################################

rule feature_counts:
    input:
        sam=expand("results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
                   sample=config["Samples"]),
        annotation=config["Reference_Genome"]["Annotation"],
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        "results/LCMSeq/gene_quanti/counts/all.featureCounts",
        "results/LCMSeq/gene_quanti/counts/all.featureCounts.summary",
        "results/LCMSeq/gene_quanti/counts/all.featureCounts.jcounts"
    threads: 6
    params:
        extra=(
            " -O --fracOverlap 0.2"
            " -t exon -g gene_id"
            " -s 1" # fr-secondstrand
            " -p" # Pair-End
        )
    log:
        config["Log_Dir"] + "/featurecounts/all.log"
    conda:
        "../envs/featurecounts.yaml"
    script:
        "../wrappers/featurecounts.py"

rule calculate_TPM:
    input:
        bam="results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
        gtf=config["Reference_Genome"]["Annotation"]
    output:
        out="results/LCMSeq/gene_quanti/tpm/{sample}.out",
        ent="results/LCMSeq/gene_quanti/tpm/{sample}.ent",
        uni="results/LCMSeq/gene_quanti/tpm/{sample}.uni"
    log: config["Log_Dir"] + "/TPMCalculator/{sample}.log"
    conda:
        "../envs/tpmcalculator.yaml"
    script:
        "../wrappers/tpmcalculator.py"
