{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--cell-order"
            },
            "type": "File",
            "id": "Cell_Order"
        },
        {
            "inputBinding": {
                "prefix": "--gene-list"
            },
            "type": "File",
            "id": "Gene_List"
        },
        {
            "inputBinding": {
                "prefix": "--sparse-data-table"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Sparse_Data_Table"
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
                "glob": "*.csv.gz"
            },
            "type": "File",
            "id": "Data_Tables"
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
        "mist_sparse_to_dense.py"
    ],
    "class": "CommandLineTool",
    "id": "SparsetoDense",
    "cwlVersion": "v1.0"
}