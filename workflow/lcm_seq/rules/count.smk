###############################################################################
# TRANSCRIPTOME QUATIFICATION
###############################################################################

rule gene_quanti:
    input:
        bam="results/align/alignments/{sample}.sorted_byname.bam",
        gtf="data/annotations/Triticum_aestivum.IWGSC.43.gtf"
    output:
        "results/gene_quanti/gene_counts_per_lib/{sample}.tsv"
    shell:
        " htseq-count -r name -i gene_id -s yes "
        " -f bam {input.bam} {input.gtf} > {output} "

rule gene_cal_tmp:
    input:
        bam="results/align/alignments/{sample}.sorted_byname.bam",
        gtf="data/annotations/Triticum_aestivum.IWGSC.43.gtf"
    output:
        out="results/gene_quanti/gene_tpm_per_lib/{sample}.out",
        ent="results/gene_quanti/gene_tpm_per_lib/{sample}.ent",
        uni="results/gene_quanti/gene_tpm_per_lib/{sample}.uni"
    run:
        import os
        import glob
        pwd = os.getcwd()
        cwd = os.path.dirname(output.out)
        gtf_file = os.path.join(pwd, input.gtf)
        bam_file = os.path.join(pwd, input.bam)
        shell(
            " cd {cwd} && TPMCalculator -v -p "
            " -g {gtf_file} -b {bam_file} -k gene_id "
        )
        all_output = glob.glob(os.path.join(cwd, "%s*" % wildcards.sample))
        for out_file in all_output:
            os.rename(out_file, out_file.replace(".sorted_byname_genes", ""))

rule diff_sample:
    '''
    Generate sample information table for DESeq2
    '''
    output:
        "results/gene_quanti/differential_expression/{contrast}/sampleTable.tsv"
    run:
        sampleTable = ["{0}\t{0}.tsv\t{1}\n" 
            .format(sample, config["SAMPLE_INFO"][sample]["condition"]) 
            for sample in config["DE_ANALYSIS"][wildcards.contrast]["sample"]]
        ref_cond = config["DE_ANALYSIS"][wildcards.contrast]["reference"]
        with open(output[0], "w") as sampleFile:
            sampleFile.writelines(sampleTable)

rule diff_gene:
    # Dependencies:
    #   - r-optparse
    #   - bioconductor-deseq2
    #   - bioconductor-apeglm
    input: lambda wildcards: [ \
        "results/gene_quanti/gene_counts_per_lib/{}.tsv".format(sample) \
            for sample in config["DE_ANALYSIS"][wildcards.contrast]["sample"]],
        sampleTable = "results/gene_quanti/differential_expression/{contrast}/sampleTable.tsv"
    output:
        lfcResults="results/gene_quanti/differential_expression/{contrast}/results.tsv",
        lfcShrink="results/gene_quanti/differential_expression/{contrast}/results_shrunken.tsv",
        normCounts="results/gene_quanti/differential_expression/{contrast}/norm_counts.tsv",
        mergedResults="results/gene_quanti/differential_expression/{contrast}/merged_results.tsv",
        mergedShrink="results/gene_quanti/differential_expression/{contrast}/merged_results_shrunken.tsv",
        rdata="results/gene_quanti/differential_expression/{contrast}/debug.RData"
    params:
        htseqdir="results/gene_quanti/gene_counts_per_lib",
        inputType = "htseq",
        threshold = 10,
        reference = lambda wildcards: config["DE_ANALYSIS"][wildcards.contrast]["reference"]
    script:
        "src/scripts/deseq2.R"

rule viz_sample_distance:
    input: "results/gene_quanti/differential_expression/{contrast}/debug.RData"
    output:
        variance_plot="results/gene_quanti/differential_expression/{contrast}/plots/variance_plot.png",
        pca_plot="results/gene_quanti/differential_expression/{contrast}/plots/pca_plot.png"
    script:
        "src/scripts/viz_sample_distance.R"

rule pull_go_data:
    '''
    The table is directly pulled from BioMart, and contains all genes.
    '''
    output: "results/gene_quanti/differential_expression/GO_DATA.tsv"
    params:
        species="Triticum_aestivum"
    script:
        "src/scripts/pull_go_data.R"


rule anno_gene:
    '''
    Pull annotation from BioMart.
    '''
    # Dependencies:
    #   - bioconductor-biomart
    #   - r-dplyr
    input: "results/gene_quanti/differential_expression/{contrast}/merged_results.tsv"
    output: "results/gene_quanti/differential_expression/{contrast}/merged_results_with_anno.tsv"
    params:
        filterCol="ensembl_gene_id",
        species="Triticum_aestivum"
    script: 
        "src/scripts/anno_gene.R"

rule go_enrich:
    # Dependencies:
    #   - bioconductor-clusterProfiler
    #   - 
    '''
    Go enrichment by using clusterProfiler.
    '''
    input: 
        DEGs="results/gene_quanti/differential_expression/{contrast}/merged_results_with_anno.tsv",
        GOData="results/gene_quanti/differential_expression/GO_DATA.tsv"
    output: 
        png="results/gene_quanti/differential_expression/{contrast}/go_enrichment.png",
        tsv="results/gene_quanti/differential_expression/{contrast}/go_enrichment.tsv"
    params:
        filterCol="ensembl_gene_id",
        species="Triticum_aestivum"
    script:
        "src/scripts/go_enrich.R"