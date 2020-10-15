{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--basic-algo-only"
            },
            "type": [
                "null",
                "boolean"
            ],
            "id": "Basic_Algo_Only"
        },
        {
            "inputBinding": {
                "prefix": "--exact-cell-count"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "Exact_Cell_Count"
        },
        {
            "inputBinding": {
                "prefix": "--full-gene-list"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Full_Genes"
        },
        {
            "inputBinding": {
                "prefix": "--gene-status",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Gene_Status_List"
        },
        {
            "inputBinding": {
                "prefix": "--mol-annot",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Molecule_Annotation_List"
        },
        {
            "inputBinding": {
                "prefix": "--putative-cell-call"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "Putative_Cell_Call"
        },
        {
            "inputBinding": {
                "prefix": "--seq-stats"
            },
            "type": "File",
            "id": "Seq_Metrics"
        },
        {
            "inputBinding": {
                "prefix": "--tag-names",
                "itemSeparator": ","
            },
            "type": [
                "null",
                {
                    "items": "string",
                    "type": "array"
                }
            ],
            "id": "Tag_Names"
        }
    ],
    "hints": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        }
    ],
    "outputs": [
        {
            "outputBinding": {
                "glob": "metrics-files.tar.gz"
            },
            "type": "File",
            "id": "Annot_Files"
        },
        {
            "outputBinding": {
                "glob": "Cell_Label_Filtering/*.png"
            },
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "id": "Cell_Label_Filter"
        },
        {
            "outputBinding": {
                "glob": "cell_order.json"
            },
            "type": "File",
            "id": "Cell_Order"
        },
        {
            "outputBinding": {
                "glob": "*_Expression_Data.st.gz"
            },
            "type": "File",
            "id": "Expression_Data"
        },
        {
            "outputBinding": {
                "glob": "*_Expression_Data_Unfiltered.st.gz"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Expression_Data_Unfiltered"
        },
        {
            "outputBinding": {
                "glob": "gene_list.json"
            },
            "type": "File",
            "id": "Gene_List"
        },
        {
            "outputBinding": {
                "glob": "Annotations/*_Annotation_Molecule.csv.gz"
            },
            "type": "File",
            "id": "Molecular_Annotation"
        },
        {
            "outputBinding": {
                "glob": "Cell_Label_Filtering/*_Putative_Cells_Origin.csv"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Putative_Cells_Origin"
        },
        {
            "outputBinding": {
                "glob": "*PerCell_Sparse.csv.gz"
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Sparse_Data_Tables"
        },
        {
            "outputBinding": {
                "glob": "*RSEC*PerCell_Unfiltered_Sparse.csv.gz"
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Sparse_Data_Tables_Unfiltered"
        },
        {
            "outputBinding": {
                "glob": "Annotations/*_Annotation_Molecule_Trueno.csv"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Tag_Annotation"
        },
        {
            "outputBinding": {
                "glob": "Trueno/*_Calls.csv"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Tag_Calls"
        },
        {
            "outputBinding": {
                "glob": "Trueno/*"
            },
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "id": "Trueno_out"
        },
        {
            "outputBinding": {
                "glob": "Annotations/*_UMI_Adjusted_CellLabel_Stats.csv"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "UMI_Adjusted_CellLabel_Stats"
        },
        {
            "outputBinding": {
                "glob": "Annotations/*_UMI_Adjusted_Stats.csv"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "UMI_Adjusted_Stats"
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
        "mist_get_datatables.py"
    ],
    "class": "CommandLineTool",
    "id": "GetDataTable",
    "cwlVersion": "v1.0"
}