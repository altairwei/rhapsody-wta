#include <iostream>
#include <set>
#include <vector>

#include "BamTools/api/BamReader.h"
#include "BamTools/api/BamWriter.h"

#include "utils.h"

using namespace std;
using namespace BamTools;

struct bucket
{
    size_t count = 0;
    set<string> genes;
};


int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "BAM file is required." << endl;
        return 1;
    }

    string inputFilename = argv[1];

    BamReader reader;
    if ( !reader.Open(inputFilename) ) {
        cerr << "Could not open input BAM files." << endl;
        return 1;
    }

    // look up index files for all BAM files
    reader.LocateIndex();

    // make sure index data is available
    /*
    if (!reader.HasIndex()) {
        cerr << "ERROR: could not load index data for all input BAM "
                "file... Aborting."
             << endl;
        reader.Close();
        return 1;
    }
    */

    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    cerr << "Start to loop BAM file..." << endl;

    int n_depth = 100;
    vector<bucket> buckets(n_depth);
    BamAlignment al;
    size_t i = 0;
    while ( reader.GetNextAlignmentCore(al) ) {
        string tag_value;
        if (al.GetTag("XF", tag_value)) {
            int idx = getRandomInt(0, n_depth);
            buckets[idx].count += 1;
            buckets[idx].genes.insert(tag_value);
        }

        cerr << "Processed " << ++i << " alignments" << "\r" << flush;
    }
    // close the reader & writer
    reader.Close();
}
