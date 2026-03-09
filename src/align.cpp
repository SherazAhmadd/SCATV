#include "aligner.h"
#include <algorithm>

void circularviral_aligner::align_sequences(string fastq_file, string output_file) {
    
    ifstream input_file(fastq_file); 
    ofstream out(output_file); 

    if (!input_file.is_open()) {
        cout << "Error: The fastq reads failed to open " << fastq_file << endl;
        exit(1);
    }
    if (!out.is_open()) {
        cout << "Error: Could not create output file " << output_file << endl;
        exit(1);
    }

    string line;
    int line_count = 0;
    string read_header = "";
    string read_seq = "";


    out << "@HD\tVN:1.6\tSO:unsorted" << endl;
    out << "@SQ\tSN:VIRAL_CIRCULAR\tLN:" << original_length << endl;
    out << "@PG\tID:SCATV\tPN:Spaced-seed_Circular_Aligner_Targeted_Virus\tVN:1.0" << endl;

    while (getline(input_file, line)) {
        line_count++;
        int line_type = line_count % 4;

        if (line_type == 1) {
            read_header = line; 
        } else if (line_type == 2) {
            read_seq = line;    
            
            int mask_len = space_mask.length();
            unordered_map<int, int> position_votes; 
            
            // Spaced Seed Indexing algo
            for (int i = 0; i <= (int)read_seq.length() - mask_len; i++) {
                string current_window = read_seq.substr(i, mask_len);
                string masked_kmer = "";
                for (int j = 0; j < mask_len; j++) {
                    masked_kmer += (space_mask[j] == '1') ? current_window[j] : '*';
                }
                
                if (spaced_index.find(masked_kmer) != spaced_index.end()) {
                    for (int matching_pos : spaced_index[masked_kmer]) {
                        int estimated_start = matching_pos - i; 
                        if (estimated_start >= 0) position_votes[estimated_start]++;
                    }
                }
            }
            
            int best_anchor = -1, max_votes = 0;
            for (auto const& [pos, votes] : position_votes) {
                if (votes > max_votes) { max_votes = votes; best_anchor = pos; }
            }

            if (best_anchor != -1) {
                // Define window with 20base pair buffer
                int window_start = max(0, best_anchor - 20);
                int window_length = min((int)read_seq.length() + 40, (int)double_reference.length() - window_start);
                string ref_window = double_reference.substr(window_start, window_length);
                
                // Smith-Waterman Alignment
                int rows = read_seq.length() + 1;
                int cols = ref_window.length() + 1;
                vector<vector<int>> score_matrix(rows, vector<int>(cols, 0));
                vector<vector<int>> traceback(rows, vector<int>(cols, 0)); 
                
                int max_sw_score = 0, max_r = 0, max_c = 0;

                for (int r = 1; r < rows; r++) {
                    for (int c = 1; c < cols; c++) {
                        int score = (read_seq[r-1] == ref_window[c-1]) ? Match_score : Mismatch_score;
                        int diag = score_matrix[r-1][c-1] + score;
                        int up = score_matrix[r-1][c] + Gap_score;
                        int left = score_matrix[r][c-1] + Gap_score;
                        
                        int best = max({0, diag, up, left});
                        score_matrix[r][c] = best;
                        if (best == 0) traceback[r][c] = 0;
                        else if (best == diag) traceback[r][c] = 1;
                        else if (best == up) traceback[r][c] = 2;
                        else traceback[r][c] = 3;

                        if (best > max_sw_score) { max_sw_score = best; max_r = r; max_c = c; }
                    }
                }

                // CIGAR Traceback
                string raw_cigar = "";
                int curr_r = max_r, curr_c = max_c;
                while (curr_r > 0 && curr_c > 0 && traceback[curr_r][curr_c] != 0) {
                    if (traceback[curr_r][curr_c] == 1) { raw_cigar += 'M'; curr_r--; curr_c--; }
                    else if (traceback[curr_r][curr_c] == 2) { raw_cigar += 'I'; curr_r--; }
                    else { raw_cigar += 'D'; curr_c--; }
                }
                reverse(raw_cigar.begin(), raw_cigar.end());
                
                string final_cigar = "";
                if (!raw_cigar.empty()) {
                    int count = 1;
                    for(int i = 1; i <= raw_cigar.length(); i++) {
                        if(i < raw_cigar.length() && raw_cigar[i] == raw_cigar[i-1]) count++;
                        else { final_cigar += to_string(count) + raw_cigar[i-1]; count = 1; }
                    }
                }

  
                int absolute_pos = window_start + curr_c + 1; 
                if (absolute_pos > original_length) absolute_pos -= original_length;
                // writing the results
    
                out << read_header.substr(1) << "\t0\tVIRAL_CIRCULAR\t" 
                    << absolute_pos << "\t255\t" << final_cigar << "\t*\t0\t0\t" 
                    << read_seq << "\t*" << endl;
            }
        } 
    }

    input_file.close();
    out.close();
    cout << "[SCATV] Alignment finished. Results saved to " << output_file << endl;

}
