## Feature `linkMerge` of cwl has broken

Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52

Error log:

```log
INFO:toil.worker:---TOIL WORKER OUTPUT LOG---
INFO:toil:Running Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52.
Traceback (most recent call last):
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/worker.py", line 366, in workerScript
    job._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore, defer=defer)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/job.py", line 1392, in _runner
    returnValues = self._run(jobGraph, fileStore)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/job.py", line 1329, in _run
    return self.run(fileStore)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 566, in run
    resolved_cwljob = resolve_indirect(self.cwljob)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 210, in resolve_indirect
    res = _resolve_indirect_inner(inner)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 179, in _resolve_indirect_inner
    result[key] = value.resolve()
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 126, in resolve
    source = promise[1][promise[0]]
TypeError: tuple indices must be integers or slices, not str
ERROR:toil.worker:Exiting the worker because of a failed job on host tanglab
WARNING:toil.jobGraph:Due to failure we are reducing the remaining retry count of job 'file:///home/testuser/src/test-folder/tests/test.cwl#test-link-merge/BundleLogs/c4e1067b-8800-4141-a2cb-d4de50f92a55' cat kind-file_home_testuser_src_test-folder_tests_test.cwl_test-link-merge_BundleLogs_c4e1067b-8800-4141-a2cb-d4de50f92a55/instance7q_s67gx with ID kind-file_home_testuser_src_test-folder_tests_test.cwl_test-link-merge_BundleLogs_c4e1067b-8800-4141-a2cb-d4de50f92a55/instance7q_s67gx to 0
```


## How to reproduce the problem

Here is the content of cwl document file `test.cwl`:

```yaml
class: Workflow
cwlVersion: v1.0
id: test-link-merge
inputs:
  - id: echo1_log
    type: string
  - id: echo2_log
    type: string
outputs:
  - id: log
    outputSource:
      - BundleLogs/logs
    type: File
requirements:
  - class: MultipleInputFeatureRequirement
steps:
  - id: echo1
    in:
      - id: echo1_op
        source: echo1_log
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand: echo
      inputs:
        - id: echo1_op
          type: string
          inputBinding:
            position: 1
      outputs:
        - id: output
          type: stdout
      stdout: echo1.log
  - id: echo2
    in:
      - id: echo2_op
        source: echo2_log
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand: echo
      inputs:
        - id: echo2_op
          type: string
          inputBinding:
            position: 1
      outputs:
        - id: output
          type: stdout
      stdout: echo2.log
  - id: BundleLogs
    in:
      - id: log_files
        linkMerge: merge_flattened
        source:
          - echo1/output
          - echo2/output
    out:
      - id: logs
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand: cat
      inputs:
        - id: log_files
          type: 'File[]'
          inputBinding:
            position: 1
      outputs:
        - id: logs
          type: stdout
      stdout: logs.txt
```

Here is the content of inputs object file `test.yml`:

```yaml
echo1_log: Hello World 1
echo2_log: Hello World 2
```

### Use toil-cwl-runner:

```shell
toil-cwl-runner --writeLogs . --logFile cwltoil.log test.cwl test.yml
```

Error log:

```log
INFO:toil.worker:---TOIL WORKER OUTPUT LOG---
INFO:toil:Running Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52.
Traceback (most recent call last):
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/worker.py", line 366, in workerScript
    job._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore, defer=defer)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/job.py", line 1392, in _runner
    returnValues = self._run(jobGraph, fileStore)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/job.py", line 1329, in _run
    return self.run(fileStore)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 566, in run
    resolved_cwljob = resolve_indirect(self.cwljob)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 210, in resolve_indirect
    res = _resolve_indirect_inner(inner)
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 179, in _resolve_indirect_inner
    result[key] = value.resolve()
  File "/home/testuser/usr/miniconda3/envs/wta/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 126, in resolve
    source = promise[1][promise[0]]
TypeError: tuple indices must be integers or slices, not str
ERROR:toil.worker:Exiting the worker because of a failed job on host tanglab
WARNING:toil.jobGraph:Due to failure we are reducing the remaining retry count of job 'file:///home/testuser/src/test-folder/tests/test.cwl#test-link-merge/BundleLogs/c4e1067b-8800-4141-a2cb-d4de50f92a55' cat kind-file_home_testuser_src_test-folder_tests_test.cwl_test-link-merge_BundleLogs_c4e1067b-8800-4141-a2cb-d4de50f92a55/instance7q_s67gx with ID kind-file_home_testuser_src_test-folder_tests_test.cwl_test-link-merge_BundleLogs_c4e1067b-8800-4141-a2cb-d4de50f92a55/instance7q_s67gx to 0
```

### Use cwltool:

```shell
cwltool test.cwl test.yml
```

The output:

```log
INFO /home/testuser/usr/miniconda3/envs/wta/bin/cwltool 1.0.20190906054215
INFO Resolved 'test.cwl' to 'file:///home/testuser/src/test-folder/tests/test.cwl'
INFO [workflow ] start
INFO [workflow ] starting step echo2
INFO [step echo2] start
INFO [job echo2] /tmp/se6vhf20$ echo \
    'Hello World 2' > /tmp/se6vhf20/echo2.log
INFO [job echo2] completed success
INFO [step echo2] completed success
INFO [workflow ] starting step echo1
INFO [step echo1] start
INFO [job echo1] /tmp/u2uimnbd$ echo \
    'Hello World 1' > /tmp/u2uimnbd/echo1.log
INFO [job echo1] completed success
INFO [step echo1] completed success
INFO [workflow ] starting step BundleLogs
INFO [step BundleLogs] start
INFO [job BundleLogs] /tmp/fl895pdf$ cat \
    /tmp/tmpv92lp2p3/stga2fda6e1-aa18-4efb-a1c4-e797a960859b/echo1.log \
    /tmp/tmpv92lp2p3/stgc5075679-9ee0-4e63-a78f-7ce48922afb4/echo2.log > /tmp/fl895pdf/logs.txt
INFO [job BundleLogs] completed success
INFO [step BundleLogs] completed success
INFO [workflow ] completed success
{
    "log": {
        "location": "file:///home/testuser/src/test-folder/tests/logs.txt",
        "basename": "logs.txt",
        "class": "File",
        "checksum": "sha1$8a32893e6f95776abbd184219bea56f23dac800e",
        "size": 28,
        "path": "/home/testuser/src/test-folder/tests/logs.txt"
    }
}
INFO Final process status is success
```