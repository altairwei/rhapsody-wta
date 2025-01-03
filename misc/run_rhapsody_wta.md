# Run BD Rhapsody WTA pipeline

## 准备用户空间 Docker 运行时

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
  --writeLogs logs \
  --logFile cwltoil.log \
  --logLevel INFO \
  --retryCount 2 \
  --maxLogFileSize 20000000000 \
  --stats \
  rhapsody-wta-yaml.cwl template_wta.yml \
&
```

使用 cwltool ：

```shell
mkdir -p tmp/docker_tmp
export TMPDIR=$(pwd)/tmp/docker_tmp
cwltool --parallel --outdir results rhapsody-wta-yaml.cwl template_wta.yml
```

如果使用 NFS 文件系统：

> Using a shared filesystem: Toil currently only supports a tempdir set to a local, non-shared directory.

也就是说，在 NFS 文件系统上不能设置 `--workDir` 。

下面是在 sippe 服务器上执行的命令：

```shell
toil-cwl-runner \
  --user-space-docker-cmd=udocker \
  --jobStore file:results/rhapsody-wta-job-store \
  --outdir results \
  --writeLogs logs \
  --logFile cwltoil.log \
  --logLevel INFO \
  --retryCount 2 \
  --maxLogFileSize 20000000000 \
  --stats \
  rhapsody-wta-yaml.cwl template_wta.yml
```

使用 cwltool :

```shell
mkdir -p tmp/docker_tmp
export TMPDIR=$(pwd)/tmp/docker_tmp
cwltool --user-space-docker-cmd=udocker --parallel --outdir results rhapsody-wta-yaml.cwl template_wta.yml
```

值得注意的是，udocker v1.1.4 只支持 python2 ，所以请安装开发版的 udocker 。但是开发版本的 udocker 有缺陷，cwltool 使用起来会报错。

### 原始测序数据文件名约定

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

### 使用额外的临时文件夹确保储存空间足够

docker 会产生很多临时文件，因此需要确保 /tmp 目录有足够空间。

也可以改变临时文件位置。

```shell
mkdir docker_tmp
export TMPDIR=`pwd`/docker_tmp
```

### 参考基因组索引使用的 STAR 版本要与 docker 镜像中的保持一致

`AnnotateR2` 这一步发生问题：

```shell
Running the following: `['STAR', '--runThreadN', '6', '--genomeDir', '.', '--readFilesIn', '/var/lib/cwl/stgbe721f1e-b127-47e4-bf8a-53493828320e/SRR6354915_1_R2_.fastq.gz', '--outSAMunmapped Within', '--outFilterScoreMinOverLread', '0', '--outFilterMatchNminOverLread', '0', '--outFilterMultimapScoreRange', '0', '--seedSearchStartLmax', '50', '--readFilesCommand', 'gunzip', '-c', '--outFileNamePrefix', 'SRR6354915_1.', '--quantMode', 'TranscriptomeSAM', '--quantTranscriptomeBan', 'Singleend', '--outSAMorder', 'PairedKeepInputOrder']`
Feb 24 08:08:52 ..... started STAR run
Feb 24 08:08:53 ..... loading genome

EXITING: FATAL INPUT ERROR: unrecoginzed parameter name "genomeFileSizes" in input "genomeParameters.txt"
SOLUTION: use correct parameter name (check the manual)
```

Docker 镜像中使用的 STAR 版本为 2.5.2b ，可能是因为基因组索引版本不兼容。

### 基因组索引文件压缩包的约定

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

### 使用 toil 的 jobstore 来实现 CWL 流程的可重入

中间结果没有被缓存下来，每次运行失败只能重新来过？

CWL 的参考实现 cwltool 文档中根本没有提可重入性的问题，`--cachedir` 也存在 BUG。

比如:

```shell
cwl-runner --tmpdir-prefix tmp/ --cachedir cache/ --outdir results/ workflow.cwl job.yml
```

而且 rhapsody-wta 流程又涉及到大量的 Docker ，要不要保留 docker 容器也是个未知的问题。

Toil 具有 jobstore 可以实现工作流程的重入。

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

~~我觉得应该是删除出错的 job ，然后再重启 jobstore~~

### 文件名太长

```
OSError: [Errno 36] File name too long
```

将整个 Workflow 的 id 改得短一点。

### BundleLogs 错误

```log
INFO:toil.worker:---TOIL WORKER OUTPUT LOG---
INFO:toil:Running Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52.
Traceback (most recent call last):
  File "/home/altairwei/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/worker.py", line 366, in workerScript
    job._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore, defer=defer)
  File "/home/altairwei/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/job.py", line 1392, in _runner
    returnValues = self._run(jobGraph, fileStore)
  File "/home/altairwei/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/job.py", line 1329, in _run
    return self.run(fileStore)
  File "/home/altairwei/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 566, in run
    resolved_cwljob = resolve_indirect(self.cwljob)
  File "/home/altairwei/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 210, in resolve_indirect
    res = _resolve_indirect_inner(inner)
  File "/home/altairwei/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 179, in _resolve_indirect_inner
    result[key] = value.resolve()
  File "/home/altairwei/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 126, in resolve
    source = promise[1][promise[0]]
TypeError: tuple indices must be integers or slices, not str
ERROR:toil.worker:Exiting the worker because of a failed job on host tanglab
WARNING:toil.jobGraph:Due to failure we are reducing the remaining retry count of job 'file:///home/altairwei/src/rhapsody-wta/rhapsody-wta-yaml.cwl#rhapsody-wta/BundleLogs/1b6dd7a2-def6-4c8c-a41c-80d8715e18a1' kind-file_home_altairwei_src_rhapsody-wta_rhapsody-wta-yaml.cwl_rhapsody-wta_BundleLogs_1b6dd7a2-def6-4c8c-a41c-80d8715e18a1/instancenjt87_px with ID kind-file_home_altairwei_src_rhapsody-wta_rhapsody-wta-yaml.cwl_rhapsody-wta_BundleLogs_1b6dd7a2-def6-4c8c-a41c-80d8715e18a1/instancenjt87_px to 0
```

`self.sources` 就等于 `to_merge` 对象。所以 `promise` 对象相当于下面这个列表中的单个元素：

```python
to_merge = [(shortname(s), promises[s].rv()) for s in aslist(source_obj)]
```

`shortname(s)` 是路径字符串的最后一个组分；而 `source_obj` 就是 CWL 中 `source` 下面的 `array<string>` 。

那么 `source = promise[1][promise[0]]` 中 `promise[1]` 是获取 Workflow 中 step 的运行结果，`promise[0]` 是 output 的标识符那么整个语句的涵义就是获取运行结果中的某个输出。

那么最有可能的原因是 `promise[1]` 不是 step 的输出，而是一个元组！也就是说，很可能是 toil 的 bug !

问题出在这：

```python
    def run(self, file_store):
        resolved_cwljob = resolve_indirect(self.cwljob)
        metadata = {}
        if isinstance(resolved_cwljob, tuple):
            cwljob = resolved_cwljob[0]
            metadata = resolved_cwljob[1]
        else:
            cwljob = resolved_cwljob
        fill_in_defaults(
            self.cwltool.tool['inputs'], cwljob,
            self.runtime_context.make_fs_access(
                self.runtime_context.basedir or ""))
        realjob = CWLJob(self.cwltool, cwljob, self.runtime_context)
        self.addChild(realjob)
        return realjob.rv(), metadata
```

多返回了一个 `metadata`！

见相关的问题：

- https://github.com/DataBiosphere/toil/issues/2846
- https://github.com/DataBiosphere/toil/pull/2845

~~是不是 `Logs` 于之间建立的 `logs` 文件夹重名了？~~

有可能时 BundleLogs 依赖的某个source没有了其实不存在，以至于报错？我怀疑不能直接引用 `XXX/output`

出问题的步骤：

```yaml
  - id: BundleLogs
    in:
      - id: log_files
        linkMerge: merge_flattened
        source:
          - AnnotateReads/output
          - AnnotateR1/output
          - AnnotateR2/output
          - CheckReference/output
          - GetDataTable/output
          - Metrics/output
          - AddtoBam/output
          - AnnotateMolecules/output
          - QualityFilter/output
          - CheckFastqs/log
          - SplitAndSubsample/log
          - MergeBAM/log
          - Sparse_to_Dense_Datatable/output
          - Sparse_to_Dense_Datatable_Unfiltered/output
          - IndexBAM/log
    out:
      - id: logs_dir
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: |-
        ${
          /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */
          function uuid() {
            var uuid = "", i, random;
            for (i = 0; i < 32; i++) {
              random = Math.random() * 16 | 0;
              if (i == 8 || i == 12 || i == 16 || i == 20) {
                uuid += "-";
              }
              uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);
            }
            return uuid;
          }
          var listing = [];
          for (var i = 0; i < inputs.log_files.length; i++) {
            var log_file = inputs.log_files[i];
            log_file.basename = uuid() + "-" + log_file.basename;
            listing.push(log_file);
          }
          return ({
            logs_dir: {
              class: "Directory",
              basename: "Logs",
              listing: listing
            }
          });
        }
      hints: []
      inputs:
        - id: log_files
          type: 'File[]'
      outputs:
        - id: logs_dir
          type: Directory
      requirements:
        - class: InlineJavascriptRequirement
        - class: MultipleInputFeatureRequirement
      
    requirements: []
```

在 run.py 中运行 cwl ，然后使用 pdb 来 debug