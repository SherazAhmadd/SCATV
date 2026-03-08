#include "aligner.h"

void print_usage() {
    cout << "++++++++ SCATV: Spaced-seed Circular Aligner ++++++++" << endl;
    cout << "Usage: ./SCATV -r <ref.fna> -q <reads.fq> -o <results.sam>" << endl;
    cout << "Options:" << endl;
    cout << "  -r    Path to the doubled reference FASTA file" << endl;
    cout << "  -q    Path to the simulated FASTQ reads" << endl;
    cout << "  -o    Path to the output SAM file" << endl;
    cout << "_______________________________________________" << endl;
}

int main(int argc, char* argv[]) {

    if (argc < 7) {
        print_usage();
        return 1;
    }

    string ref_file = "";
    string query_file = "";
    string output_file = "";

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        
        if (arg == "-r" && i + 1 < argc) {
            ref_file = argv[++i];
        } else if (arg == "-q" && i + 1 < argc) {
            query_file = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            output_file = argv[++i];
        } else {
            cout << "Unknown or incomplete argument: " << arg << endl;
            print_usage();
            return 1;
        }
    }

    if (ref_file == "" || query_file == "" || output_file == "") {
        cout << "All flags (-r, -q, and -o) are required." << endl;
        print_usage();
        return 1;
    }

    cout << "Starting Alignment..." << endl;
    cout << "Reference File: " << ref_file << endl;
    cout << "Query File:     " << query_file << endl;
    cout << "Output File:    " << output_file << endl;

    circularviral_aligner my_aligner;
    
    my_aligner.read_fasta(ref_file);
    my_aligner.build_index();
    
    my_aligner.align_sequences(query_file, output_file);

    cout << "Alignment finished. Results saved to " << output_file << endl;

    return 0;
}