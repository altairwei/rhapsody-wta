{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--extra-seqs"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Extra_Seqs"
        },
        {
            "inputBinding": {
                "prefix": "--index",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Index"
        },
        {
            "inputBinding": {
                "prefix": "--R2"
            },
            "type": "File",
            "id": "R2"
        }
    ],
    "requirements": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        },
        {
            "envDef": [
                {
                    "envName": "CORES_ALLOCATED_PER_CWL_PROCESS",
                    "envValue": "$(String(runtime.cores))"
                }
            ],
            "class": "EnvVarRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "outputs": [
        {
            "outputBinding": {
                "glob": "*Annotation_R2.csv.gz"
            },
            "type": "File",
            "id": "Annot_R2"
        },
        {
            "outputBinding": {
                "glob": "*-annot.gtf"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "GTF"
        },
        {
            "outputBinding": {
                "glob": "*mapping_R2.BAM"
            },
            "type": "File",
            "id": "R2_Bam"
        },
        {
            "outputBinding": {
                "glob": "*_picard_quality_metrics.csv.gz"
            },
            "type": "File",
            "id": "R2_Quality_Metrics"
        },
        {
            "outputBinding": {
                "glob": "*.log"
            },
            "type": "File",
            "id": "output"
        }
    ],
    "baseCommand": [
        "mist_annotate_R2.py"
    ],
    "class": "CommandLineTool",
    "id": "AnnotateR2",
    "cwlVersion": "v1.0"
}