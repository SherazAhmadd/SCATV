#include "aligner.h"

void circularviral_aligner::build_index() {
    
    int ref_len = double_reference.length();
    int mask_len = space_mask.length();
    
    if (ref_len == 0) {
        cout << "Error: Reference sequence is empty. Cannot build index." << endl;
        exit(1);
    }

    cout << "Building spaced-seed index with mask: " << space_mask << endl;

    // Slide the window across the entire doubled reference
    for (int i = 0; i <= ref_len - mask_len; i++) {
        
        string current_window = double_reference.substr(i, mask_len);
        string masked_kmer = "";
        
        // Apply the 110110111 mask to the current window
        for (int j = 0; j < mask_len; j++) {
            if (space_mask[j] == '1') {
                masked_kmer += current_window[j];
            } else {
                masked_kmer += '*'; 
            }
        }
        
        // Add the starting coordinate to our hash map
        spaced_index[masked_kmer].push_back(i);
    }

    cout << "Index successfully built. Total unique spaced seeds: " << spaced_index.size() << endl;
}