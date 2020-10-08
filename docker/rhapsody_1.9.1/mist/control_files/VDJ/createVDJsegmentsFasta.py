from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def loadCleanFasta(fastafile):
    # Downstream tools can choke on certain characters in fasta headers

    recordList = list(SeqIO.parse(fastafile, 'fasta'))

    for record in recordList:
        record.id = record.id.replace('(', '_') \
                             .replace(')', '_') \
                             .replace('[', '_') \
                             .replace(']', '_')

    return recordList

# Want to balance the need for long enough sequences that 200+bp reads will align successfully, with the desire for speed
# More sequences in the reference will slow down the alignment
# V segments are long enough to be aligned against on their own
# D, J, and C region segments are not always long enough - so concat these in all possible combinations

# Note: IGKappa locus has a 'distal' inverted V segments that can result in inverted J segments

# Note: Our IGLC primer will also amplify IGLL5 mRNA. (Uses the same constant region exon) We add into the reference these sequences as a mRNA target.



species_list = ["HomoSapiens", "MusMusculus"]

cFragments = {"HomoSapiens": { "TRAC": "CAGATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTG",
                               "TRBC": "AGGACCTGAANAANGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATC",
                               "TRGC": "NATAAACAACTTGATGCAGATGTTTCCC",
                               "TRDC": "GAAGTCAGCCTCATACCAAACCATCCGTTTTTGTCATGAAAAATGGAACAAATGTCGCTTGTCTGGTGAAGGAATTCTACCCCAAGGATAT",
                               "IGKC": "GAACTGTGGCTGCACCATCTGTCTTCATCTTCCCGCCATCTGA",
                               "IGLC": "GTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCNCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGT",
                               "IGHM": "GGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGT",
                               "IGHD": "CACCCACCAAGGCTCCGGATGTGTTCCCCATCATATCAGGGTGCAGACA",
                               "IGHA": "CATCCCCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTCNNCAGCACCCNNCNAGATGGGAACGTGGTCNTCGCNTGCCTGGTCCAGGGCTTCTTCCCCCAGGAGCCACTCAGTGTGACCTGGAGCGAAAG",
                               "IGHG": "CNTCCACCAAGGGCCCATCGGTCTTCCCCCTGGCNCCCTNCTCCANGAGCACCTCNGNGNGCACAGCNGCCCTGGGCTGCCTGGTCAAGGACTACTT",
                               "IGHE": "CCTCCACACAGAGCCCATCCGTCTTCCCCTTGAC"
                            },
               "MusMusculus": { "TRAC": "NACATCCAGAACCCAGAACCT",
                                "TRBC": "AGGATCTGAGAAATGTGACTCCACCCAAGGTCTCCTTGTTTGAGCCATCAAAAGCAGAGATTG",
                                "TRGC": "ACAAANGNNNTGANNCAGACNTTTCNCCCAAGCCTACTATT",
                                "TRDC": "GCTTGTCTGGTGAAAGATTTCTAC",
                                "IGKC": "CNACTGTATCCATCTTCCCACCATCCAGTGAGCAGTTAACATCT",
                                "IGLC": "GNCAGCCCAAGTCNNCNCCNNCNNTCACCNTGTTTCCACCTTCCNCTGANGAGCTC",
                                "IGHM": "AGAGTCAGTCCTTCCCAAATGTC",
                                "IGHD": "GTGATAAAAAGGAACCTGACATGTTCCTCCTCTCAGAGTGCAAAGCCCCAGAGGAAAATGAAAAGATAAACCTGGGCTGTTTAGTAATTGGAAGTCAG",
                                "IGHA": "AGNCTGCNAGANANCCCACCATCTACCCACTGACA",
                                "IGHG": "CNANAACNACACCCCCATCNGTCTATCCNNTGGNCCCTG",
                                "IGHE": "TCTATCAGGAACCCTCAGCTCTA"
                            }
             }

v_seg_dict = {"TR": {"A":"AlphaT", "B":"BetaT", "G":"GammaT", "D":"DeltaT"},
            "IG": {"K":"KappaB", "L":"LambdaB", "H":"HeavyB"}
            }

light_seg_dict = {"TR": {"A":"AlphaT", "G":"GammaT"},
                    "IG": {"K":"KappaB", "L":"LambdaB"}
                    }

heavy_seg_dict = {"TR": {"B":"BetaT", "D":"DeltaT"},
                "IG": {"H":"HeavyB"}
                }

heavyTypes = ["M", "D", "A", "G", "E"]


for species in species_list:

    AllSeqs = {"TR": {}, "IG": {}}
    SeqRecords = {"TR": [], "IG": []}

    # V genes
    for majorType in v_seg_dict:

        for chain in v_seg_dict[majorType]:

            vFasta = loadCleanFasta(species + "_" + majorType + chain + "V.fasta")
            for vgene in vFasta:
                revCompSeq = vgene.seq.upper().reverse_complement()
                if str(revCompSeq) in AllSeqs[majorType]:
                    print("Skipping duplicateSeq:" + vgene.id)
                else:
                    SeqRecords[majorType].append(SeqRecord(revCompSeq, id=vgene.id + "|rc_" + v_seg_dict[majorType][chain] + "VDJ", description="", name=""))
                    AllSeqs[majorType][str(revCompSeq)] = 1


    # DJC gene combinations - Heavy
    for majorType in heavy_seg_dict:

        for chain in heavy_seg_dict[majorType]:

            dFasta = loadCleanFasta(species + "_" + majorType + chain + "D.fasta")
            jFasta = loadCleanFasta(species + "_" + majorType + chain + "J.fasta")

            if majorType == "IG" and chain == 'H':
                cFragmentNames = [ majorType + chain + hType for hType in heavyTypes ]
            else:
                cFragmentNames = [ majorType + chain + "C" ]

            cFragmentDict = {cGene:cFragments[species][cGene] for cGene in cFragmentNames}

            for dgene in dFasta:
                for jgene in jFasta:
                    for cFragment in cFragmentDict:
                        concat_seq = dgene.seq + jgene.seq + cFragmentDict[cFragment]
                        revComp_concat = concat_seq.upper().reverse_complement()
                        if str(revComp_concat) in AllSeqs[majorType]:
                            print("Skipping duplicateSeq:" + dgene.id + "|" + jgene.id)
                        else:
                            SeqRecords[majorType].append(SeqRecord(revComp_concat, id=dgene.id + "|" + jgene.id + "|" + cFragment + "|rc_" + heavy_seg_dict[majorType][chain] + "VDJ", description="", name=""))
                            AllSeqs[majorType][str(revComp_concat)] = 1

    # DJC gene combinations - light
    for majorType in light_seg_dict:

        for chain in light_seg_dict[majorType]:

            jFasta = loadCleanFasta(species + "_" + majorType + chain + "J.fasta")

            cFragmentName = majorType + chain + "C"
            cFragmentDict = {cFragmentName:cFragments[species][cFragmentName]}

            for jgene in jFasta:
                for cFragment in cFragmentDict:
                    concat_seq = dgene.seq + jgene.seq + cFragmentDict[cFragment]
                    revComp_concat = concat_seq.upper().reverse_complement()
                    if str(revComp_concat) in AllSeqs[majorType]:
                        print("Skipping duplicateSeq:" + jgene.id)
                    else:
                        SeqRecords[majorType].append(SeqRecord(revComp_concat, id=jgene.id + "|" + cFragment + "|rc_" + light_seg_dict[majorType][chain] + "VDJ", description="", name=""))
                        AllSeqs[majorType][str(revComp_concat)] = 1


    # if species == "HomoSapiens":
    #     ENST00000526893rc = Seq("GTCAATGAGGATATTTATTGGGGTTTCATGAGTGCAGGGAGAAGGGCTGGATGACTTGGGATGGGGAGAGAGACCCCTCCCCTGGGATCCTGCAGCTCCAGGCTCCCGTGGGTGGGGTTAGAGTTGGGAACCTATGAACATTCTGTAGGGGCCACTGTCTTCTCCACGGTGCTCCCTTCATGCGTGACCTGGCAGCTGTAGCTTCTGTGGGACTTCCACTGCTCGGGCGTCAGGCTCAGGTAGCTGCTGGCCGCGTACTTGTTGTTGCTCTGTTTGGAGGGTTTGGTGGTCTCCACTCCCGCCTTGACGGGGCTGCCATCTGCCTTCCAGGCCACTGTCACAGCTCCCGGGTAGAAGTCACTGATCAGACACACTAGTGTGGCCTTGTTGGCTTGGAGCTCCTCAGAGGAGGGCGGGAACAGAGTGACAGTGGGGTTGGCCTTGGGCTGACCTAGGACGGTGACCTTGGTCCCAGTTCCGAAGACATAACACAGTGACTGAGGCTCAGACCAAAACCCCCGGGGCCAGCACCTGGGGTCTGCTCTCTGGGGGCTGGGCTGGAGCAGGAGCCTGCCCCACAGGCTCCGCAGGCTGGATCGGCTGCTTCCAACTGAGGCTCCAGGGTCTGGGTCCCCGCTTTGCGGTGCAACCATTGGGCGCAGCAGGCCATGGGCGACCATGGCCAGACCCAGCAGCAGCAGGGGCCAGCGCTGCCTGG")
    #     SeqRecords["IG"].append(SeqRecord(ENST00000526893rc, id="IGLL5|ENST00000526893.5|rcIGLC_capture", description="", name=""))
    #     ENST00000531372rc = Seq("GTCAATGAGGATATTTATTGGGGTTTCATGAGTGCAGGGAGAAGGGCTGGATGACTTGGGATGGGGAGAGAGACCCCTCCCCTGGGATCCTGCAGCTCCAGGCTCCCGTGGGTGGGGTTAGAGTTGGGAACCTATGAACATTCTGTAGGGGCCACTGTCTTCTCCACGGTGCTCCCTTCATGCGTGACCTGGCAGCTGTAGCTTCTGTGGGACTTCCACTGCTCGGGCGTCAGGCTCAGGTAGCTGCTGGCCGCGTACTTGTTGTTGCTCTGTTTGGAGGGTTTGGTGGTCTCCACTCCCGCCTTGACGGGGCTGCCATCTGCCTTCCAGGCCACTGTCACAGCTCCCGGGTAGAAGTCACTGATCAGACACACTAGTGTGGCCTTGTTGGCTTGGAGCTCCTCAGAGGAGGGCGGGAACAGAGTGACAGTGGGGTTGGCCTTGGGCTGACCTGCCCCACAGGCTCCGCAGGCTGGATCGGCTGCTTCCAACTGAGGCTCCAGGGTCTGGGTCCCCGCTTTGCGGTGCAACCATTGGGCGCAGCAGGCCATGGGCGACCATGGCCAGACCCAGCAGCAGCAGGGGCCAGCGCTGCCTGGGACCAG")
    #     SeqRecords["IG"].append(SeqRecord(ENST00000531372rc, id="IGLL5|ENST00000531372.1|rcIGLC_capture", description="", name=""))

    for majorType in SeqRecords:
        SeqIO.write(SeqRecords[majorType], species + "_" + majorType + "_VDJsegments" + ".fasta", "fasta-2line")


