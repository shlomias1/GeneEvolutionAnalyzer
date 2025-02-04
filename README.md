# Genome Analysis Tool

## Overview
This Python program is designed for analyzing genomic data retrieved from GenBank. It extracts gene information, compares gene sets, calculates synonymous codon positions, aligns sequences, and evaluates evolutionary selection pressure using dN/dS ratios.

## Features
- Fetch genomic data from GenBank using accession numbers
- Extract gene names and compare gene sets
- Count the total and protein-coding genes
- Identify synonymous codon positions
- Remove stop codons and pad sequences for alignment
- Align sequences using a codon-based approach
- Calculate dN/dS ratios to infer selection pressure

## Dependencies
This program requires the following Python libraries:
- `Biopython`
- `pandas`

To install the dependencies, run:
```sh
pip install biopython pandas
```

## Code Structure
The program is structured into multiple classes:

### 1. `GenBankHandler`
Handles fetching genomic data from GenBank.

### 2. `GenomeAnalysis`
- Retrieves gene information.
- Compares gene sets across genomes.

### 3. `CodonAnalysis`
- Calculates synonymous codon positions.

### 4. `SequenceProcessor`
- Removes stop codons.
- Pads sequences for codon alignment.

### 5. `SequenceAligner`
- Aligns sequences using the `PairwiseAligner` from `Biopython`.

### 6. `EvolutionaryAnalysis`
- Extracts gene sequences from records.
- Calculates dN/dS ratios for evolutionary selection analysis.

## Usage
1. Set the email for GenBank queries.
2. Define a list of accession numbers for genomic sequences.
3. The program will fetch data, process sequences, and display results in tabular format.

## Output
- A table of synonymous codon positions.
- Gene statistics for each genome.
- A comparative analysis of shared and unique genes.
- A dN/dS evolutionary comparison table.

## Example Execution
Run the program with:
```sh
python genome_analysis.py
```

Example output:
```
===== Synonymous Positions =====
Codon  | Synonymous Positions
----------------------------
ATA    | 2
ATC    | 2
...

===== Gene Statistics =====
GenBank Data for NC_045512.2:
Total number of genes: 23
Number of protein-coding genes: 12
...

===== Evolutionary Comparison =====
Gene       |  dN   |  dS   | dN/dS |  Selection type
--------------------------------------------------
S          | 2.520 | 1.671 | 1.508 | Positive selection
ORF7b      | 0.010 | 0.039 | 0.248 | Negative selection
...
```

## License
This project is licensed under the MIT License.
