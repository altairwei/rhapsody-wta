Firstly, we need to setup environments:

```shell
source env.sh
```

How to run workflow:

```shell
run_toil.sh workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```

Restart last workflow without losing any progress:

```shell
run_toil.sh -r workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```

Start the same workflow from scratch:

```shell
run_toil.sh -c workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```

Run toil leader within a named screen session:

```shell
run_toil.sh -s workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```