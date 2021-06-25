# Generating Genome Indexes

## 使用 STAR 构建人类基因组的索引

```shell
nohup \
STAR \
    --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir data/reference_indexes \
    --genomeFastaFiles data/reference_sequences/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile data/annotations/Homo_sapiens.GRCh38.99.gtf \
&
```

当索引构建好后就不在需要原始参考基因组序列了。

## 将索引和参考基因组打包

> For WTA assays, the reference genome is a compressed tarball that contains the **STAR index files** for the species of the cells used in the BD WTA experiment.

基因组索引文件需要放置到一个文件夹下，然后使用打包成 tarball。这个文件夹会被用来当做参考基因组的名称。

```shell
# 首先创建压缩包
nohup \
tar --transform 's#^\.#GRCh38#' \
    -cvf data/reference_sequences/Homo_sapiens.GRCh38.tar.gz -C data/reference_indexes . \
&
```

可以使用 `--use-compress-program=pigz` 选项来并行压缩，不过似乎格式不兼容？

主要涉及到如下流程：

```yaml
class: CommandLineTool
cwlVersion: v1.0
baseCommand:
- mist_check_references.py
inputs:
- id: Reference
    type: 'File[]'
    inputBinding:
    position: 0
    prefix: '--reference'
    itemSeparator: ','
    shellQuote: false
- id: Label_Version
    type: int?
    inputBinding:
    position: 0
    prefix: '--label-version'
    shellQuote: false
- id: AbSeq_Reference
    type: 'File[]?'
    inputBinding:
    position: 0
    prefix: '--abseq-reference'
    shellQuote: false
- id: Sample_Tags_Version
    type: string?
    inputBinding:
    position: 0
    prefix: '--sample-tags-version'
    shellQuote: false
- id: Putative_Cell_Call
    type: int?
    inputBinding:
    position: 0
    prefix: '--putative-cell-call'
    shellQuote: false
outputs:
- id: Index
    type: 'File[]'
    outputBinding:
    glob: '*-annot.*'
    outputEval: |
        ${
            if (self.length == 1) { // Targeted
                return self;
            } else if (self.length == 0){ // WTA
                return inputs.Reference;
            }
        }
- id: Extra_Seqs
    type: File?
    outputBinding:
    glob: combined_extra_seq.fasta
- id: Reference_Panel_Names
    type: File
    outputBinding:
    glob: reference_panel_names.json
- id: output
    type: File
    outputBinding:
    glob: '*.log'
- id: Full_Genes
    type: File?
    outputBinding:
    glob: full-gene-list.json
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
        dockerPull: 'images.sbgenomics.com/aberno/rhapsody:1.8'
    - class: InlineJavascriptRequirement
```