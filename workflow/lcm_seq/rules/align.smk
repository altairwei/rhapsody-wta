###############################################################################
# READS MAPPING
###############################################################################

rule build_refseq_index:
    input: config["GENOME_REF_SEQ"]["Triticum_aestivum"]
    output: expand("results/align/index/{species}.{num}.ht2l", \
        species="Triticum_aestivum", num=list(range(1,9)))
    threads: 16
    params:
        #FIXME: bmax is not obvious, this should be replaced by human-readable unit.
        bmax="306940781",
        index_basename="results/align/index/Triticum_aestivum"
    shell:
        "hisat2-build -p {threads} --bmax {params.bmax} --dcv 1024 -a" 
        " -f --large-index {input} {params.index_basename} "

rule hisat2_mapping:
    input: 
        index_files=expand("results/align/index/{species}.{num}.ht2l", \
            species="Triticum_aestivum", num=list(range(1,9))),
        forward="results/quality_control/trimmed_reads/{sample}_1P.fq.gz",
        reverse="results/quality_control/trimmed_reads/{sample}_2P.fq.gz",
        #FIXME: if unpaired reads were mapped to reference, ht-seq should not use `-r name` option.
        #unpaired_forward="results/quality_control/trimmed_reads/{sample}_1U.fq.gz",
        #unpaired_reverse="results/quality_control/trimmed_reads/{sample}_2U.fq.gz"
    params:
        index_basename="results/align/index/Triticum_aestivum"
    output:
        summary="results/align/alignments/{sample}_summary.log",
        matched=temp("results/align/alignments/{sample}.sam"),
        unmatched="results/align/unaligned_reads/{sample}_unmatched.fq.gz"
    threads: 8
    shell:
        "hisat2 -q -p {threads} -x {params.index_basename} "
        " --rna-strandness FR "
        " -1 {input.forward} -2 {input.reverse} " 
        #" -U {input.unpaired_forward},{input.unpaired_reverse} "
        " -S {output.matched} --un-gz {output.unmatched} "
        " --summary-file {output.summary} "

rule sort_bam_by_pos:
    input: 
        "results/align/alignments/{sample}.sam",
    output:
        "results/align/alignments/{sample}.sorted_bypos.bam",
    params:
        prefix="results/align/alignments/{sample}.sorted_bypos"
    threads: 4
    shell:
        "samtools sort -@ {threads} -T {params.prefix} -O bam -o {output} {input}"

rule sort_bam_by_name:
    input: 
        "results/align/alignments/{sample}.sorted_bypos.bam",
    output:
        "results/align/alignments/{sample}.sorted_byname.bam",
    params:
        prefix="results/align/alignments/{sample}.sorted_byname"
    threads: 4
    shell:
        "samtools sort -n -@ {threads} -T {params.prefix} -O bam -o {output} {input}"