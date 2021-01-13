#include <iostream>
#include <set>
#include <vector>
#include <cstdio>

#include "BamTools/api/BamReader.h"
#include "BamTools/api/BamWriter.h"
#include "CLI/CLI.hpp"

#include "utils.h"

using namespace std;
using namespace BamTools;

struct bucket
{
    size_t read_count = 0;
    set<string> genes;
};


int main(int argc, char *argv[]) {
    string inputFilename;
    bool verbose;

    CLI::App app{"Saturation"};
    app.add_option("file", inputFilename, "SAM file")->required();
    app.add_flag("-v,--verbose", verbose, "Show progress");
    CLI11_PARSE(app, argc, argv);

    BamReader reader;
    if ( !reader.Open(inputFilename) ) {
        cerr << "Could not open input BAM files." << endl;
        return 1;
    }

    // look up index files for all BAM files
    reader.LocateIndex();

    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    int n_depth = 100;
    vector<bucket> buckets(n_depth);
    BamAlignment al;
    size_t i = 0;
    string tag_value;
    while ( reader.GetNextAlignment(al) ) {
        // Select a bucket for each reads.
        int idx = getRandomRange(0, n_depth);
        buckets[idx].read_count += 1;

        if (al.GetTag<string>("XF", tag_value))
            buckets[idx].genes.insert(tag_value);

        if (verbose) {
            if (++i % 1000000 == 0)
                fprintf(stderr, "Processed %zu alignments\r", i);
        }


    }

    reader.Close();

    // Write header of CSV ouput
    cout << "depths" << "," << "detected_genes" << endl;

    size_t depth = 0;
    set<string> detected_genes;
    for (auto &bucket : buckets) {
        depth += bucket.read_count;
        detected_genes.insert(bucket.genes.begin(), bucket.genes.end());
        //cout << bucket.read_count << "," << bucket.genes.size() << endl;
        cout << depth << "," << detected_genes.size() << endl;
    }
}
