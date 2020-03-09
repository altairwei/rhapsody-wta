# Run BD Rhapsody WTA pipeline

## 准备用户控件 Docker 运行时

```shell
conda install singularity -n wta
```

## 试运行

这个工作流程无法使用 singularity ，不过可以尝试下 udocker 或者 dx-docker

```shell
mkdir -p tmp/docker_tmp
export TMPDIR=$(pwd)/tmp/docker_tmp
nohup \
toil-cwl-runner \
  --jobStore file:results/rhapsody-wta-job-store \
  --workDir tmp \
  --outdir results \
  --writeLogs results/logs \
  --logFile cwltoil.log \
  --logLevel INFO \
  --retryCount 2 \
  --maxLogFileSize 20000000000 \
  --stats \
  rhapsody-wta-yaml.cwl template_wta.yml \
&
```

如果使用 NFS 文件系统：

> Using a shared filesystem: Toil currently only supports a tempdir set to a local, non-shared directory.

也就是说，再 NFS 文件系统上不能设置 `--workDir` 。

下面是再 sippe 服务器上执行的命令：

```shell
toil-cwl-runner \
  --user-space-docker-cmd=udocker \
  --jobStore file:results/rhapsody-wta-job-store \
  --outdir results \
  --writeLogs results/logs \
  --logFile cwltoil.log \
  --logLevel INFO \
  --retryCount 2 \
  --maxLogFileSize 20000000000 \
  --stats \
  rhapsody-wta-yaml.cwl template_wta.yml
```

### 文件名不符合规则

文件名：

```yaml
Reads:

 - class: File
   location: "data/reads/SRR6354915_R1_.fastq.gz"

 - class: File
   location: "data/reads/SRR6354915_R2_.fastq.gz"
```

官方文档说明：

- An underscore on each side of R1 or R2 (_R1_ and _R2_).
- The <sample> name should be the same for R1 and R2.
- Convert uncompressed files to .gz format.
- **Basename of file should end with 001.**

Example:

    <sample>_S1_L001_R1_001.fastq.gz
    <sample>_S1_L001_R2_001.fastq.gz

工作流程 `PairReadFiles`这一步用来匹配名称的正则表达式：

```javascript
// This is the illumina basespace regex. More sophisticated file handling is needed for NovaSeq
// example: <SAMPLE>[<SAMPLE NUMBER>][<LANE>]_R<READ FLAG>_001.fastq.gz
/^(.*?)(_S[0-9]*)?(_L[0-9]*)?(_R[1|2])_001(-[0-9]*)?\.(.*?)$/
```

> ### Naming Convention
> FASTQ files are named with the sample name and the sample number, which is a numeric assignment based on the order that the sample is listed in the sample sheet. For example: Data\Intensities\BaseCalls\SampleName_S1_L001_R1_001.fastq.gz
>
> - SampleName — The sample name provided in the sample sheet. If a sample name is not provided, the file name includes the sample ID, which is a required field in the sample sheet and must be unique.
> - S1 — The sample number based on the order that samples are listed in the sample sheet starting with 1. In this example, S1 indicates that this sample is the first sample listed in the sample sheet. **NOTE** Reads that cannot be assigned to any sample are written to a FASTQ file for sample number 0, and excluded from downstream analysis.
>
> - L001 — The lane number.
> - R1 — The read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.
> - 001 — The last segment is always 001.
>
> FASTQ files that do not follow this naming convention cannot be imported into BaseSpace.

### 参考基因组索引问题

`AnnotateR2` 这一步发生问题：

```shell
Running the following: `['STAR', '--runThreadN', '6', '--genomeDir', '.', '--readFilesIn', '/var/lib/cwl/stgbe721f1e-b127-47e4-bf8a-53493828320e/SRR6354915_1_R2_.fastq.gz', '--outSAMunmapped Within', '--outFilterScoreMinOverLread', '0', '--outFilterMatchNminOverLread', '0', '--outFilterMultimapScoreRange', '0', '--seedSearchStartLmax', '50', '--readFilesCommand', 'gunzip', '-c', '--outFileNamePrefix', 'SRR6354915_1.', '--quantMode', 'TranscriptomeSAM', '--quantTranscriptomeBan', 'Singleend', '--outSAMorder', 'PairedKeepInputOrder']`
Feb 24 08:08:52 ..... started STAR run
Feb 24 08:08:53 ..... loading genome

EXITING: FATAL INPUT ERROR: unrecoginzed parameter name "genomeFileSizes" in input "genomeParameters.txt"
SOLUTION: use correct parameter name (check the manual)
```

Docker 镜像中使用的 STAR 版本为 2.5.2b ，可能是因为基因组索引版本不兼容。

### 磁盘空间不够

docker 会产生很多临时文件，因此需要确保 /tmp 目录有足够空间。

也可以改变临时文件位置。

```shell
mkdir docker_tmp
export TMPDIR=`pwd`/docker_tmp
```

### 基因组索引文件错误

重点在于 `OSError: [Errno 22] Invalid argument: '.'` 这段，我怀疑是打包基因组索引时引入的路径前缀 `.`。

```shell
Runtime.totalMemory()=536870912
...done
Cleaning up...
Traceback (most recent call last):
  File "src/mist/lib/MistLogger.py", line 32, in mist.lib.MistLogger.log_death.log_death_inner
  File "src/mist/apps/AnnotateR2.py", line 118, in mist.apps.AnnotateR2.annotate_r2
  File "src/mist/apps/utils.py", line 449, in mist.apps.utils.cleanup
  File "/opt/conda/lib/python3.6/shutil.py", line 484, in rmtree
    onerror(os.rmdir, path, sys.exc_info())
  File "/opt/conda/lib/python3.6/shutil.py", line 482, in rmtree
    os.rmdir(path)
OSError: [Errno 22] Invalid argument: '.'

INFO [job AnnotateR2] Max memory used: 29731MiB
ERROR [job AnnotateR2] Job error:
("Error collecting output for parameter 'Annot_R2':\nrhapsody-wta-json.cwl:296:25: Did not find output file with glob pattern: '['*Annotation_R2.csv.gz']'", {})
WARNING [job AnnotateR2] completed permanentFail
WARNING [step AnnotateR2] completed permanentFail
```

不过在去除掉 `.` 前缀后，又出现问题。重点在这 `'--genomeDir', 'chrLength.txt'` ，这说明前缀文件夹是必要的！

```shell
Mapping read 2 via STAR...
Running the following: `['STAR', '--runThreadN', '6', '--genomeDir', 'chrLength.txt', '--readFilesIn', '/var/lib/cwl/stg00e26ac7-c14d-40dd-b503-299b0116f532/SRR6354915_1_R2_.fastq.gz', '--outSAMunmapped Within', '--outFilterScoreMinOverLread', '0', '--outFilterMatchNminOverLread', '0', '--outFilterMultimapScoreRange', '0', '--seedSearchStartLmax', '50', '--readFilesCommand', 'gunzip', '-c', '--outFileNamePrefix', 'SRR6354915_1.', '--quantMode', 'TranscriptomeSAM', '--quantTranscriptomeBan', 'Singleend', '--outSAMorder', 'PairedKeepInputOrder']`
Feb 25 02:52:17 ..... started STAR run
Feb 25 02:52:17 ..... loading genome

EXITING because of FATAL ERROR: could not open genome file chrLength.txt/genomeParameters.txt
SOLUTION: check that the path to genome files, specified in --genomeDir is correct and the files are present, and have user read permsissions
```

查看一下源代码：

```python
        # untar the STAR index - find archive independently of name
        with tarfile.open(star_idx_tarball) as tar:
            star_idx = tar.getnames()[0]
            tar.extractall()

        ref_name = star_idx.split('.')[0]
```

看起来基因组索引需要放在一个文件夹下面，而这个文件夹会被用来当做参考基因组名字。

### 内存空间不足

`GetDataTable` 需要 64G 内存空间。

```shell
ERROR Exception on step 'GetDataTable'
ERROR [step GetDataTable] Cannot make job: Requested at least 64000 ram but only 62559 available
```

对内存资源的要求写在了流程里：

```yaml
    requirements:
      - class: ResourceRequirement
        ramMin: 64000
```

直接先改成 32G 试试。

`Metrics` 也需求很大的内存空间：

```shell
INFO [workflow ] starting step Metrics
INFO [step Metrics] start
ERROR Exception on step 'Metrics'
ERROR [step Metrics] Cannot make job: Requested at least 96000 ram but only 62410 available
```

在改变掉内存需求后，总算成功了！

### 为什么CWL流程无法重入？

中间结果没有被缓存下来，每次运行失败只能重新来过？

CWL 的参考实现 cwltool 文档中根本没有提可重入性的问题，`--cachedir` 也存在 BUG。

比如:

```shell
cwl-runner --tmpdir-prefix tmp/ --cachedir cache/ --outdir results/ workflow.cwl job.yml
```

而且 rhapsody-wta 流程又涉及到大量的 Docker ，要不要保留 docker 容器也是个未知的问题。

### 工作流程错误

下面是错误输出：

```shell
("Error collecting output for parameter 'Cell_Order':\nrhapsody-wta-yaml.cwl:1533:13: Did not find output file with glob pattern: '['random_stdout_89ff98c4-e6a8-44f2-9619-4cb8e11955c6']'", {})
```

出错的流程：

```yaml
  - id: Sparse_to_Dense_File
    in:
      - id: GDT_cell_order
        source: GetDataTable/Cell_Order
    out:
      - id: Cell_Order
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand:
        - cat
      inputs:
        - id: GDT_cell_order
          type: File?
          inputBinding:
            position: 1
            shellQuote: false
      outputs:
        - id: Cell_Order
          type: File
          outputBinding:
            glob: random_stdout_89ff98c4-e6a8-44f2-9619-4cb8e11955c6
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
      stdout: cell_order.json
    requirements: []
```

按照最新的 cwl 文档，应该改成：

```yaml
      outputs:
        - id: Cell_Order
          type: stdout
```

我要如何让 toil 依据新的 cwl 文件重启运行呢？

我觉得应该是删除出错的 job ，然后再重启 jobstore