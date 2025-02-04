# Evolutionary Genomic Analysis Tool

## Overview
This Python program analyzes evolutionary relationships between genomes using dN/dS ratios. It fetches genomic data from GenBank, extracts gene sequences, aligns them, and calculates synonymous and non-synonymous mutation rates. The program also generates a comparative bar plot for dN and dS values.

## Features
- Fetches GenBank data for specified accession numbers.
- Extracts genome length, gene counts, and protein-coding genes.
- Identifies common and unique genes between genomes.
- Calculates synonymous positions for codons.
- Aligns gene sequences and removes stop codons.
- Computes dN/dS ratios to determine selection pressure.
- Visualizes results with a bar chart displaying dN and dS values.

## Dependencies
Ensure you have the following Python libraries installed:
```bash
pip install biopython matplotlib pandas numpy
```

## Usage
### Running the Program
Execute the script with:
```bash
python script.py
```

### Expected Output
- **Synonymous Positions**: A table displaying synonymous counts for each codon.
- **Gene Statistics**: Total and protein-coding gene counts for each genome.
- **Gene Comparison**: Lists common and unique genes.
- **Evolutionary Analysis**: A table with dN, dS, dN/dS ratios, and selection type.
- **Graph**: A bar chart comparing dN and dS values for each gene.

### Example Output
```
===== Synonymous Positions =====
Codon  | Synonymous Positions
ATA    | 2
ATC    | 2
...
===== Gene Statistics =====
GenBank Data for NC_045512.2:
Total number of genes: 23
Number of protein-coding genes: 12
...
===== Evolutionary Comparison =====
Gene  |  dN   |  dS   | dN/dS | Selection type
--------------------------------------------------
S      | 2.520 | 1.671 | 1.508 | Positive selection
...
```

## Code Structure
### Classes
#### 1. `GenBankHandler`
Handles fetching genomic data from GenBank.
#### 2. `GenomeAnalysis`
Processes genome features including gene extraction and comparison.
#### 3. `CodonAnalysis`
Calculates synonymous positions for codons.
#### 4. `SequenceProcessor`
Removes stop codons and pads sequences.
#### 5. `SequenceAligner`
Performs global sequence alignment.
#### 6. `EvolutionaryAnalysis`
Calculates dN/dS ratios and plots comparative graphs.

### Functionality
- `fetch_data()`: Retrieves genome data from GenBank.
- `calculate_dn_ds()`: Computes dN, dS, and dN/dS ratios.
- `plot_dn_ds()`: Plots a bar chart of dN and dS values with labels.

## Graphical Output
- The generated graph compares dN and dS values for each gene.
- Each bar displays its respective value above it.

## Future Improvements
- Support additional evolutionary models.
- Allow user-defined input for accession numbers.
- Save results to a file for further analysis.

## Author
Developed by Shlomi Assayag.
