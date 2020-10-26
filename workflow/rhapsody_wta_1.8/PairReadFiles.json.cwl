{
    "inputs": [
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Reads"
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "outputs": [
        {
            "type": {
                "items": {
                    "fields": [
                        {
                            "type": "File",
                            "name": "R1"
                        },
                        {
                            "type": "File",
                            "name": "R2"
                        },
                        {
                            "type": "int",
                            "name": "readPairID"
                        },
                        {
                            "type": "string",
                            "name": "library"
                        }
                    ],
                    "type": "record"
                },
                "type": "array"
            },
            "id": "ReadPairs"
        }
    ],
    "class": "ExpressionTool",
    "doc": "PairReadsFiles takes an array of split files and pairs them, such that an R1 file is transferred to a QualityFilter with its corresponding R2 file.\nEach file should be formatted as illumina outputs it from basespace: e.g. sample_L001_R1_001.fastq.gz. After being split, that sample file would be an array files named sample_L001_R1_001-00.fastq, sample_L001_R1_001-01.fastq, etc\n",
    "expression": "${\n  // send paired reads to the same key in readPairs\n  var readPairs = {}\n  for (var i = 0; i < inputs.Reads.length; i++) {\n    var f = inputs.Reads[i];\n\n    // This is the illumina basespace regex. More sophisticated file handling is needed for NovaSeq\n    // example: <SAMPLE>[<SAMPLE NUMBER>][<LANE>]_R<READ FLAG>_001.fastq.gz\n    var groups = f.basename.match(/^(.*?)(_S[0-9]*)?(_L[0-9]*)?(_R[1|2])_001(-[0-9]*)?\\.(.*?)$/);\n    var library = groups[1];\n    var sampleNumber = groups[2];\n    var laneNumber = groups[3];\n    var flag = groups[4];\n    var chunkID = 9999; // if there is no scatter id, use an arbitrary number\n    if (groups[5]){\n      chunkID = parseInt(groups[5].slice(1)); // slice off the '-'\n    }\n\n    // double check we have a chunk id\n    if (chunkID === undefined || chunkID === null) {\n          throw new Error(\"chunkID could not be determined!\");\n    }\n\n    // notice, we ignore the flag. This causes the paired reads to share the same key\n    var readPairID = library + sampleNumber + laneNumber + chunkID\n\n    // sort the information from the filename into an object\n    if (!readPairs[readPairID]) {\n      readPairs[readPairID] = {\n        R1: null,\n        R2: null,\n        library: library,\n        readPairID: null,\n      };\n    }\n    // add in the readPair, depending on the flag\n    if (flag === \"_R1\") {\n      readPairs[readPairID].R1 = f\n    } else if (flag === \"_R2\") {\n      readPairs[readPairID].R2 = f\n    }\n\n  }\n  // we are not interested in the keys in readPairs; flatten into an array of objects\n  var readPairsList = []\n  var i = 1;\n  for (var key in readPairs) {\n    if (readPairs.hasOwnProperty(key)) {\n      var readPair = readPairs[key]\n      readPair.readPairID = i\n      readPairsList.push(readPair)\n      i++;\n    }\n  }\n  // pass this array to the record array named \"ReadPairs\" on the CWL layer\n  return {ReadPairs: readPairsList}\n}",
    "id": "PairReadFiles",
    "cwlVersion": "v1.0"
}