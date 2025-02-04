# Gene Evolution Analyzer

## Overview
This project analyzes the evolution of genes by comparing nucleotide sequences from different GenBank records. It calculates dN/dS ratios to determine the selection pressure on genes and identifies common genes between two genomes.

## Features
- **Fetch GenBank Data**: Retrieves nucleotide sequences and annotations from GenBank.
- **Gene Comparison**: Identifies common and unique genes between two genomes.
- **Sequence Alignment**: Aligns gene sequences using pairwise codon-based alignment.
- **dN/dS Calculation**: Computes synonymous (dS) and non-synonymous (dN) substitution rates.
- **Stop Codon Removal**: Cleans sequences by removing stop codons.
- **Selection Type Analysis**: Determines whether a gene is under positive, negative, or neutral selection based on dN/dS ratio.

## Requirements
To run this project, ensure you have the following installed:
- Python 3.x
- Biopython (`pip install biopython`)

## Usage
1. Clone the repository or download the script.
2. Modify the `accession_numbers` list in `main()` to specify the genomes to compare.
3. Run the script:
   ```bash
   python GeneEvolutionAnalyzer.py
   ```
4. The script will fetch data from GenBank, process the sequences, and output a table with dN, dS, and dN/dS ratios for common genes.

## Output Example
```
Gene       |  dN   |  dS   | dN/dS |  Selection type
--------------------------------------------------
S          | 2.520 | 1.671 | 1.508 | Positive selection
ORF7b      | 0.010 | 0.039 | 0.248 | Negative selection
N          | 1.496 | 1.716 | 0.872 | Negative selection
ORF10      | 0.000 | 0.000 |  inf  | Positive selection
ORF8       | 0.000 | 0.000 |  inf  | Positive selection
```

## Functions Explained
### `fetch_genbank_data(accession_number)`
Fetches GenBank data and returns a list of `SeqRecord` objects.

### `compare_gene_sets(gene_set_1, gene_set_2, genome_1, genome_2)`
Finds common and unique genes between two genomes.

### `calculate_dn_ds_for_common_genes(common_genes, records_1, records_2)`
Aligns common gene sequences and calculates dN, dS, and dN/dS values.

### `remove_stop_codons(seq)`
Removes stop codons (`TAA`, `TAG`, `TGA`) from a given sequence.

## Notes
- The script requires an active internet connection to fetch GenBank data.
- Stop codons are removed before processing sequences.
- If dS is `0`, dN/dS is set to `inf` (infinity).

## License
This project is open-source under the MIT License.
