# Memory Consumption of Each Steps

## 对内存有要求的步骤

```yaml
  - AnnotateR2
    requirements:
      - class: ResourceRequirement
        coresMin: 6
        ramMin: 48000

  - AnnotateReads
    requirements:
      - class: ResourceRequirement
        ramMin: 32000

  - AnnotateMolecules
    requirements:
      - class: ResourceRequirement
        ramMin: 32000

  - GetDataTable
    requirements:
      - class: ResourceRequirement
        ramMin: 64000

  - Sparse_to_Dense_Datatable
    requirements:
      - class: ResourceRequirement
        ramMin: 4000

  - Sparse_to_Dense_Datatable_Unfiltered
    requirements:
      - class: ResourceRequirement
        ramMin: 4000

  - AddtoBam
    requirements:
      - class: ResourceRequirement
        outdirMin: 131072
        ramMin: 16000
        tmpdirMin: 262144

  - Metrics
    requirements:
      - class: ResourceRequirement
        outdirMin: 65536
        ramMin: 96000
        tmpdirMin: 65536
```

手动将所有超过 64G 的需求改成 48G 。

## 第一次实时运行测试

输入文件：

```shell
用量 6.8G
3.3G SRR6354915_R1_001.fastq.gz
3.5G SRR6354915_R2_001.fastq.gz
```

- AnnotateR2: 
- AnnotateReads
- AnnotateMolecules
- GetDataTable
- Sparse_to_Dense_Datatable
- Sparse_to_Dense_Datatable_Unfiltered
- AddtoBam
- Metrics

通过 `grep 'Max memory used' nohup.out` 来查看日志：

```
INFO [job CheckFastqs] Max memory used: 32MiB
INFO [job CheckReference] Max memory used: 32MiB
INFO [job SplitAndSubsample] Max memory used: 160MiB
INFO [job SplitAndSubsample_2] Max memory used: 160MiB
INFO [job QualityFilter_4] Max memory used: 89MiB
INFO [job QualityFilter] Max memory used: 102MiB
INFO [job QualityFilter_2] Max memory used: 102MiB
INFO [job QualityFilter_3] Max memory used: 102MiB
INFO [job AnnotateR1_4] Max memory used: 57MiB
INFO [job AnnotateR1_2] Max memory used: 64MiB
INFO [job AnnotateR1] Max memory used: 64MiB
INFO [job AnnotateR1_3] Max memory used: 64MiB
INFO [job AnnotateR2] Max memory used: 30181MiB
INFO [job AnnotateR2_2] Max memory used: 661MiB
INFO [job AnnotateR2_3] Max memory used: 256MiB
INFO [job AnnotateR2_4] Max memory used: 250MiB
INFO [job AnnotateReads] Max memory used: 205MiB
INFO [job AnnotateMolecules] Max memory used: 147MiB
INFO [job GetDataTable] Max memory used: 192MiB
INFO [job Sparse_to_Dense_File] Max memory used: 0MiB
INFO [job Sparse_to_Dense_Datatable] Max memory used: 64MiB
INFO [job Sparse_to_Dense_Datatable_2] Max memory used: 64MiB
INFO [job Sparse_to_Dense_Datatable_3] Max memory used: 64MiB
INFO [job Sparse_to_Dense_Datatable_4] Max memory used: 64MiB
INFO [job Sparse_to_Dense_Datatable_Unfiltered_2] Max memory used: 64MiB
INFO [job Sparse_to_Dense_Datatable_Unfiltered] Max memory used: 64MiB
INFO [job Uncompress_Expression_Matrix] Max memory used: 0MiB
INFO [job Uncompress_Datatable_2] Max memory used: 0MiB
INFO [job Uncompress_Datatable] Max memory used: 0MiB
INFO [job AddtoBam_2] Max memory used: 77MiB
INFO [job AddtoBam_3] Max memory used: 969MiB
INFO [job AddtoBam] Max memory used: 976MiB
INFO [job AddtoBam_4] Max memory used: 931MiB
INFO [job Metrics] Max memory used: 1277MiB
INFO [job MergeBAM] Max memory used: 102MiB
INFO [job IndexBAM] Max memory used: 19MiB
```

通过试跑下来，`GetDataTable` 和 `Metrics` 并没有消耗那么多内存。可能得用完整数据集来测试了。

## 第二次实时运行测试

```shell
总用量 85G
3.3G 2月  24 16:00 SRR6354915_R1_001.fastq.gz
3.5G 2月  24 16:00 SRR6354915_R2_001.fastq.gz
9.9G 2月  28 09:31 SRR6354916_R1_001.fastq.gz
 11G 2月  28 09:31 SRR6354916_R2_001.fastq.gz
2.4G 2月  28 09:32 SRR6354917_R1_001.fastq.gz
2.6G 2月  28 09:32 SRR6354917_R2_001.fastq.gz
 11G 2月  28 09:32 SRR6354918_R1_001.fastq.gz
 12G 2月  28 09:32 SRR6354918_R2_001.fastq.gz
8.1G 2月  28 09:32 SRR6354919_R1_001.fastq.gz
8.9G 2月  28 09:32 SRR6354919_R2_001.fastq.gz
6.3G 2月  28 09:32 SRR6354920_R1_001.fastq.gz
6.8G 2月  28 09:32 SRR6354920_R2_001.fastq.gz
```

## 内存需求分析

### GetDataTable

以下几个步骤可能用到大量内存：

解压缩文件：

```python
    # start building metrics archive by adding an uncompressed version of the seqmetrics file
    seq_metrics_uncompressed_fp = path.join(metrics_archive_dir, path.splitext(path.basename(seq_metrics))[0])
    with gzip.open(seq_metrics, "rt") as seq_metrics_compressed, \
            open(seq_metrics_uncompressed_fp, "wt") as seq_metrics_uncompressed:
        seq_metrics_uncompressed.write(seq_metrics_compressed.read())
```

排序：

这里指定了 18G 的内存缓冲区。

```python
def sort_molAnnot_by_cell_then_gene(molAnnot, sorted_molAnnot):
    e = dict(os.environ)
    e['LC_ALL'] = 'C'
    cmd = "sort -S 18G -t , -k 1,1 -k 3,3 -k 2,2 -o {} {}".format(sorted_molAnnot, molAnnot).split(' ')
    subprocess.check_call(cmd, env=e)
```

### Metrics

