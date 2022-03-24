###############################################################################
# READS MAPPING
###############################################################################

rule align:
    input:
        fq1="results/LCMSeq/trimmed_reads/{sample}_R1.fq.gz",
        fq2="results/LCMSeq/trimmed_reads/{sample}_R2.fq.gz",
        index=config["Reference_Genome"]["Index"],
    output:
        "results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
        "results/LCMSeq/align/alignments/{sample}/ReadsPerGene.out.tab",
    log:
        config["Log_Dir"] + "/star/{sample}.log",
    params:
        index=lambda wc, input: input.index,
        extra="--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within "
            "--quantMode GeneCounts --sjdbGTFfile {}".format(
                config["Reference_Genome"]["Annotation"]
            ),
    threads: 24
    resources:
        mem=83751862272 # 78G
    conda: "../envs/star.yaml"
    script:
        "../wrappers/star_align.py"

rule sort_bam_by_name:
    input: 
        "results/LCMSeq/align/alignments/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        "results/LCMSeq/align/alignments/{sample}/Aligned.sortedByName.out.bam",
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads: 8
    conda:
        "../envs/samtools.yaml"
    script:
        "../wrappers/samtools_sort.py"