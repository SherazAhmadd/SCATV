/*
   This is the Main file for the circular viral genome aligner project.

   This file defines the main data structures and functions for:
   - Reading reference sequences
   - Building the spaced-seed index
   - Aligning query sequences from a virus simulated FASTQ reads
*/


#ifndef ALIGNER_H
#define ALIGNER_H

#include<iostream>
#include<string>
#include<vector>
#include <unordered_map> 
#include <fstream>
using namespace std;

string get_reverse_complement(string seq);
const int Match_score = 2;
const int Mismatch_score = -1;
const int Gap_score = -2;

struct alignmentresult
{
    string aligned_reference;
    string aligned_query;
    int score;
    int final_position;

};

struct circularviral_aligner

{

    string double_reference;
    int original_length;
    string space_mask = "110110111";
    unordered_map<string, vector<int>> spaced_index;

    string read_fasta(string filename);
    void build_index();
    void align_sequences(string fastq_file);

    void align_sequences(string fastq_file, string output_file);

    string get_reverse_complement(string seq);

};

#endif