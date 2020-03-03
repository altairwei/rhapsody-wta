# BD Rhapsody Sequence Reads Exploration

## 数据来源

> Yoon, Se-Jin et al. 2019. “Reliability of Human Cortical Organoid Generation.” Nature Methods 16(1): 75–78.

仅仅从 ENA 下载了 SRP126273 中的 SRR6354915


## Reads Quality and Features

BD Rhapsody 官方手册写明了质量控制的步骤：

1. Read pairs with *low sequencing quality* are first removed.

### fastp 工具

使用 `fastp` 来评估 reads 质量。对于 R1 来说，低质量的 reads 应该直接丢弃，而不是trim头尾。

```
mkdir -p results/quality_control
mkdir -p results/quality_control/trimmed_reads
fastp \
    -h results/quality_control/quality.html \
    -j results/quality_control/quality.json \
    -i data/reads/SRR6354915_R1_.fastq.gz \
    -I data/reads/SRR6354915_R2_.fastq.gz \
    -o results/quality_control/trimmed_reads/SRR6354915_R1_trimmed.fastq.gz \
    -O results/quality_control/trimmed_reads/SRR6354915_R2_trimmed.fastq.gz

seqkit seq -s data/reads/SRR6354915_R1_.fastq.gz | fa-less
```

`fastp` 无法分别针对 R1 或者 R2 进行最小长度控制。

#### Read 1 Features

R1 的质量和碱基内容如下面两张图所示，从第二幅图可以清楚的看到

![图一](./index_files/Before_filtering_read1_quality.png)

从图一可以看到，76 位有个碱基质量跌破谷底，而从图二看起来 76 位正好是 polyT 结束的位置。

![图二](./index_files/Before_filtering_read1_base_contents.png)

从图二可以清楚的看到 9 - 20 区间是 L1 片段 `ACTGGCCTGCGA`，共 12 mer，从 29 - 40 区间是 L2 片段 `GGTAGCGGTGAC` 共 12 mer。从 57 位开始则是 polyT 。值得注意的是，这与官方文档上说的不相符，L1 是 12 mer，但 L2 应该有 13 mer 才对。

值得注意的是，BD 的 Cell Label 有着不同版本：

> Specify which version of the cell label you are using: 
> (1) for 8mer, (2) for 9mer (default), (3) for Precise targeted, (4) for Precise WTA.

我们从 `AnnotateR1.py` 代码中可以看出：

```python
def get_default_sections(label):
    """ return the default starting position of each section, including the CL, linker, UMI and polyT.
        also return the reference sequences from cellkeys, as well as the polyT cutoff """

    if label == 1:
        start_pos = [0, 8, 20, 28, 40, 48, 56, 64]
        appended_startpos = [0, 20, 40, 48, 56, 64]

        refs = [[str(ref) for ref in cellkeys.cell_key1[:96]], [str(ref) for ref in cellkeys.linker1],
                [str(ref) for ref in cellkeys.cell_key2[:96]], [str(ref) for ref in cellkeys.linker2],
                [str(ref) for ref in cellkeys.cell_key3[:96]]]

    elif label == 2:
        start_pos = [0, 9, 21, 30, 43, 52, 60, 68]
        appended_startpos = [0, 21, 43, 52, 60, 68]

        refs = [[str(ref) for ref in cellkeys.v2_cell_key1[:96]], [str(ref) for ref in cellkeys.linker1],
                [str(ref) for ref in cellkeys.v2_cell_key2[:96]], [str(ref) for ref in cellkeys.linker2],
                [str(ref) for ref in cellkeys.v2_cell_key3[:96]]]

    # CL1 + L1
    cl1_l1 = [ref + refs[1][0] for ref in refs[0]]

    # CL2 + L2
    if label == 1:
        cl2_l2 = [ref + refs[3][0] for ref in refs[2]]
    else:
        cl2_l2 = [str(ref + refs[3][0] + 'A') for ref in refs[2]]

    # CL3 alone
    cl3 = refs[4]

    appended_refs = [cl1_l1, cl2_l2, cl3]

    return start_pos, appended_startpos, refs, appended_refs
```

对于版本1来说：8mer CLS1 + 12mer L1 + 8mer CLS2 + 12mer L2 + 8mer CLS3
对于版本2来说：9mer CLS1 + 12mer L1 + 9mer CLS2 + 13mer L2 + 9mer CLS3

**也就是说，使用数据集时要注意 Cell Barcode 版本！**

#### Read 2 Features

R2 在 51 位也存在一个质量下滑。

![图三](./index_files/Before_filtering_read2_quality.png)

![图四](./index_files/Before_filtering_read2_base_contents.png)

### Trimmomatic 工具

好像也不能针对单独的 R1 或者 R2 控制最小长度。