#include "aligner.h"

string circularviral_aligner::read_fasta(string filename) {
    
    ifstream input_file(filename); 
    
    if (!input_file.is_open()) {
        cout << "The reference fasta sequence failed to open" << filename << endl;
        exit(1);
    }

    string line;
    string pure_sequence = "";

    while (getline(input_file, line)) {
     
        if (line.length() == 0) {
            continue;
        }

        if (line[0] == '>') {
            continue;
        }

        for (int i = 0; i < line.length(); i++) {
            char base = line[i];
            if (base != ' ' && base != '\n' && base != '\r') {
                pure_sequence += base;
            }
        }
    }

    input_file.close();
    
    double_reference = pure_sequence;
    original_length = double_reference.length() / 2;
    
    return double_reference;
}

string circularviral_aligner::get_reverse_complement(string seq) {
    string rc = "";
    for (int i = seq.length() - 1; i >= 0; i--) {
        char base = toupper(seq[i]);
        if (base == 'A') rc += 'T';
        else if (base == 'T') rc += 'A';
        else if (base == 'C') rc += 'G';
        else if (base == 'G') rc += 'C';
        else rc += 'N';
    }
    return rc;
}