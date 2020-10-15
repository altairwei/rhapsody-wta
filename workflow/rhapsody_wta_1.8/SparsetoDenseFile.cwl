{
    "inputs": [
        {
            "inputBinding": {
                "position": 1
            },
            "type": [
                "null",
                "File"
            ],
            "id": "GDT_cell_order"
        }
    ],
    "requirements": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        }
    ],
    "stdout": "cell_order.json",
    "outputs": [
        {
            "type": "File",
            "id": "Cell_Order",
            "outputBinding": {
                "glob": "cell_order.json"
            }
        }
    ],
    "baseCommand": "cat",
    "class": "CommandLineTool",
    "id": "SparsetoDenseFile",
    "cwlVersion": "v1.0"
}