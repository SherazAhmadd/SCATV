# SCATV: Spaced-seed Circular Aligner for Targeted Viruses
SCATV (Spaced-seed Circular Aligner for Targeted Viruses) is a high-performance C++ program designed to address the specific computational challenges of aligning sequencing reads to circular viral genomes. By utilizing a double-reference strategy and spaced-seed indexing, SCATV ensures optimal alignment across genome junctions with high sensitivity to polymorphisms and sequencing errors.

In viral genomics, many pathogens (e.g., HPV, HBV, and various bacteriophages) possess circular DNA. Standard linear alignment tools often fail to map reads that "wrap around" the junction where the genome sequence ends and begins ($L \to 1$). This leads to artificial coverage drops at the most critical regions of the viral architecture.
SCATV solves this through a three-tier algorithmic approach:
1. **Reference Doubling:** The tool conceptually treats the genome as a $2L$ linear sequence to capture junction-crossing reads seamlessly.
2. **Spaced-Seed Indexing:** Instead of contiguous k-mers, SCATV uses a non-consecutive mask (110110111). This allows the initial seeding phase to remain sensitive even if a mutation occurs at a "0" position in the mask.
3. **Dynamic Programming:** Once an anchor is found, a localized Smith-Waterman alignment is performed to resolve the optimal CIGAR string, supporting Insertions, Deletions (Indels), and Mismatches.

## Key Features
* **Circular Awareness:** Handles reads overlapping the genome junction via a doubled-reference approach.
* **Spaced-Seed Indexing:** Uses a custom mask (`110110111`) for fast, mutation-sensitive mapping.
* **Smith-Waterman Alignment:** Local alignment for high-accuracy Indel/Mismatch detection.

## Project Structure
```
SCATV/
├── CMakeLists.txt        # Build system configuration
├── README.md             # Technical documentation
├── .gitignore            # Version control filters
├── include/
│   └── aligner.h         # Class definitions and scoring constants
├── src/
│   ├── main.cpp          # CLI parsing and workflow management
│   ├── parser.cpp        # FASTA/FASTQ I/O logic
│   ├── index.cpp         # Spaced-seed hash table implementation
│   └── align.cpp         # SW-Engine and Circular math
└── example/
    ├── HPV16_circ.fna      # Sample circular viral reference
    └── query_reads_1.fq     # Sample simulated viral reads
```

## Installation and Setup
### Pre-requisites
- A C++ compiler supporting C++17 
- `cmake` version 3.10 or higher
### Build Instructions
To compile the high-performance binary, follow the out-of-source build procedure:

###  Clone the repository
```bash

git https://github.com/SherazAhmadd/SCATV-Spaced-seed-Circular-Aligner-for-Targeted-Viruses.git
```
```
cd SCATV-Spaced-seed-Circular-Aligner-for-Targeted-Viruses.git
```
### Configure and Compile

```bash

mkdir build && cd build
cmake ..
make
```

## Example Run
SCATV utilizes explicit flags for data input. Users must provide a reference file, a query read file, and a destination for the SAM results.
```
./SCATV -r ../example/HPV16_circ.fna -q ../example/query_reads_1.fq -o results.sam
```
Argument Parameters:
- -r: Path to the viral reference FASTA file.
- -q: Path to the sequencing reads in FASTQ format.
- -o: Destination path for the generated SAM file.

### Note: Scoring Parameters
The alignment engine uses a localized scoring system defined in include/aligner.h:
 - Match Reward: +2
 - Mismatch Penalty: -1
 - Gap Penalty: -2

## Author and Contact
**Mr. Rana Sheraz Ahmad**  
*Integrative Omics and Molecular Modeling Laboratory, Department of Bioinformatics and Biotechnology, Government College University Faisalabad (GCUF), Faisalabad, 38000, Pakistan*  
Email: [ranasheraz.202101902@gcuf.edu.pk](mailto:ranasheraz.202101902@gcuf.edu.pk)

