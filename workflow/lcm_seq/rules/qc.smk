###############################################################################
# RAW DATA QUALITY CONTROL
###############################################################################

rule fastqc_report:
    input: "data/reads/{sample}.fq.gz"
    output:
        "results/LCMSeq/quality_control/fastqc_report/{sample}_fastqc.zip",
        "results/LCMSeq/quality_control/fastqc_report/{sample}_fastqc.html",
        directory("results/LCMSeq/quality_control/fastqc_report/{sample}_fastqc")
    params:
        outdir="results/LCMSeq/quality_control/fastqc_report"
    shell:
        # fastqc has to wait I/O, so "--threads" flag is not useful.
        "fastqc --extract -o {params.outdir} {input}"

rule multiqc_report:
    input: 
        expand("results/LCMSeq/quality_control/fastqc_report/{sample}.R{pair}_fastqc.html", \
            sample=DATASETS, pair=[1,2])
    output:
        outdir=directory("results/LCMSeq/quality_control/multiqc_report")
    params:
        indir="results/LCMSeq/quality_control/fastqc_report"
    shell:
        "multiqc -o {output.outdir} {params.indir}"

rule trim_reads:
    input: 
        "data/reads/{sample}.R1.fq.gz",
        "data/reads/{sample}.R2.fq.gz"
    output:
        "results/LCMSeq/quality_control/trimmed_reads/{sample}_1P.fq.gz",
        "results/LCMSeq/quality_control/trimmed_reads/{sample}_1U.fq.gz",
        "results/LCMSeq/quality_control/trimmed_reads/{sample}_2P.fq.gz",
        "results/LCMSeq/quality_control/trimmed_reads/{sample}_2U.fq.gz"
    params:
        adapter=config["MISC"]["adapter"]
    threads: 4
    shell:
        "trimmomatic PE -threads {threads} {input} {output} "
        "ILLUMINACLIP:{params.adapter}:2:30:10:2:true SLIDINGWINDOW:4:20 "
        "TRAILING:20 MINLEN:20 "

rule trim_all_reads:
    input: expand("results/LCMSeq/quality_control/trimmed_reads/{sample}_{pair}{marker}.fq.gz", \
        sample=DATASETS, pair=[1,2], marker=["P", "U"])
    shell:
        "printf 'Reads library has been trimmed: %s\n' {input}"

rule fastqc_trimmed_reads:
    input: "results/LCMSeq/quality_control/trimmed_reads/{sample}.fq.gz"
    output: 
        "results/LCMSeq/quality_control/trimmed_reads_qcreport/{sample}_fastqc.zip",
        "results/LCMSeq/quality_control/trimmed_reads_qcreport/{sample}_fastqc.html",
        directory("results/LCMSeq/quality_control/trimmed_reads_qcreport/{sample}_fastqc")
    params:
        outdir="results/LCMSeq/quality_control/trimmed_reads_qcreport"
    shell:
        "fastqc --extract -o {params.outdir} {input}"

rule multiqc_trimmed_reads:
    input: 
        expand("results/LCMSeq/quality_control/trimmed_reads_qcreport/{sample}_{pair}{marker}_fastqc.html", \
            sample=DATASETS, pair=[1,2], marker=["P", "U"])
    output:
        outdir=directory("results/LCMSeq/quality_control/trimmed_reads_multiqc_report")
    params:
        indir="results/LCMSeq/quality_control/trimmed_reads_qcreport"
    shell:
        "multiqc -o {output.outdir} {params.indir}"

###############################################################################
# Alignments Quality Control
###############################################################################

rule gff3_to_gtf:
    input: "data/annotations/{species}.gff3"
    output: "data/annotations/{species}.gtf"
    shell: "gffread {input} -T -o {output}"

rule bam_qc:
    input: 
        bam = "results/LCMSeq/align/alignments/{sample}.sorted_bypos.bam",
        gtf = "data/annotations/Triticum_aestivum.IWGSC.43.gtf"
    output: 
        "results/LCMSeq/align/quality_control/bam_qc/{sample}/qualimapReport.html"
    params:
        outdir="results/LCMSeq/align/quality_control/bam_qc/{sample}"
    threads: 8
    shell: 
        "qualimap bamqc -gff {input.gtf} -ip -os -nt {threads} --java-mem-size=8G "
        "-outformat html -outdir {params.outdir} -bam {input.bam} "


rule rna_seq_qc:
    input: 
        bam = "results/LCMSeq/align/alignments/{sample}.sorted_byname.bam",
        gtf = "data/annotations/Triticum_aestivum.IWGSC.43.gtf"
    output:
        report="results/LCMSeq/align/quality_control/rna_seq_qc/{sample}/qualimapReport.html", 
        counts="results/LCMSeq/align/quality_control/rna_seq_qc/{sample}/{sample}_counts.tsv",
    params:
        outdir="results/LCMSeq/align/quality_control/rna_seq_qc/{sample}"
    threads: 8 # this tool can't set threads, but it uses multiple cores acturelly.
    shell:
        #FIXME: Should notice how to deal with multi-mapped reads in wheat.
        " mkdir -p ./.tmp && "
        " JAVA_OPTS='-Djava.io.tmpdir=./.tmp' "
        " qualimap rnaseq -s -outformat html -pe --java-mem-size=8G "
        " -bam {input.bam} -gtf {input.gtf} -outdir {params.outdir} -oc {output.counts} "
        " && printf '%s\t%s\t%s\t%s\n' {wildcards.sample} Condition_{wildcards.sample} {output.counts} '2' "
        " > {output.metadata}"

rule multi_bamqc_report:
    input:
        expand("results/LCMSeq/align/quality_control/bam_qc/{sample}/qualimapReport.html", sample=DATASETS)
    output:
        outdir=directory("results/LCMSeq/align/quality_control/bamqc_multi_report")
    params:
        indir="results/LCMSeq/align/quality_control/bam_qc"
    shell:
        "multiqc -o {output.outdir} {params.indir}"

rule multi_rnaseq_qc_report:
    input:
        expand("results/LCMSeq/align/quality_control/rna_seq_qc/{sample}/qualimapReport.html", sample=DATASETS)
    output:
        outdir=directory("results/LCMSeq/align/quality_control/rnaseq_qc_multi_report")
    params:
        indir="results/LCMSeq/align/quality_control/rna_seq_qc"
    shell:
        "multiqc -o {output.outdir} {params.indir}"

rule gen_counts_qc_samples:
    input:
        expand("results/LCMSeq/align/quality_control/rna_seq_qc/{sample}/qualimapReport.html", sample=DATASETS),
    output:
        "results/LCMSeq/align/quality_control/counts_qc/samples_data.txt"
    params:
        prefix="results/LCMSeq/align/quality_control/rna_seq_qc"
    shell:
        """
        for sample in {DATASETS}; do
            printf '%s\t%s\t%s\t%s\n' $sample ${{sample%_*}} {params.prefix}/$sample/${{sample}}_counts.tsv '2' >> {output} 
        done
        """

rule gen_multi_bam_qc_samples:
    input:
       expand("results/LCMSeq/align/quality_control/bam_qc/{sample}/qualimapReport.html", sample=DATASETS)
    output:
        "results/LCMSeq/align/quality_control/multi_bam_qc/samples_data.txt"
    params:
        prefix="results/LCMSeq/align/quality_control/bam_qc"
    shell:
        """
        for sample in {DATASETS}; do
            printf '%s\t%s\t%s\t%s\n' $sample {params.prefix}/$sample/${{sample}}_counts.tsv ${{sample%_*}} >> {output}
        done
        """

rule counts_qc:
    input: 
        "results/LCMSeq/align/quality_control/counts_qc/samples_data.txt"
    output:
        "results/LCMSeq/align/quality_control/counts_qc/GlobalReport.html",
        "results/LCMSeq/align/quality_control/counts_qc/ComparisonReport.html",
        expand("results/LCMSeq/align/quality_control/counts_qc/{sample}Report.html", sample=DATASETS)
    params:
        outdir="results/LCMSeq/align/quality_control/counts_qc"
    shell:
        # Here is a bug: qualimap of conda version, haven't cp share/scripts to bin/scripts,
        # so, qualimap cannot find Rscript.
        " mkdir -p ./.tmp && "
        " JAVA_OPTS='-Djava.io.tmpdir=./.tmp' "
        " qualimap counts -c -d {input} --java-mem-size=8G "
        " -outformat html -outdir {params.outdir} "

rule multi_bam_qc:
    input: "results/LCMSeq/align/quality_control/multi_bam_qc/samples_data.txt"
    output: 
        "results/LCMSeq/align/quality_control/multi_bam_qc/multisampleBamQcReport.html"
    params:
        outdir="results/LCMSeq/align/quality_control/multi_bam_qc"
    shell:
        " mkdir -p ./.tmp && "
        " JAVA_OPTS='-Djava.io.tmpdir=./.tmp' "
        " qualimap multi-bamqc -d {input} -outdir {params.outdir} -outformat html"
     