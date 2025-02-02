
```markdown
# Virus Gene Evolution Analyzer

## Overview
The **Virus Gene Evolution Analyzer** is a tool designed to compare viral genomes and analyze the evolutionary processes shaping their genes. This tool computes the **dN/dS ratio** for common protein-coding genes shared between two viral genomes (in this case, SARS-CoV-2 strains). The **dN/dS ratio** is a key indicator of the evolutionary pressures acting on genes, classifying them into three categories:

- **Positive selection** (dN/dS > 1): Indicates that mutations that alter the protein sequence are favored, possibly suggesting adaptive evolution.
- **Negative selection** (dN/dS < 1): Indicates that mutations that alter the protein sequence are detrimental and selected against.
- **Neutral selection** (dN/dS ≈ 1): Indicates that mutations occur but do not have a significant effect on the organism's fitness.

## Features
- **Fetch Genomic Data**: Retrieves viral genomic data from GenBank using accession numbers.
- **Gene Comparison**: Identifies common and unique genes between two genomes.
- **dN/dS Calculation**: Computes the **dN (non-synonymous substitutions)** and **dS (synonymous substitutions)** for common genes between the genomes and calculates the **dN/dS ratio**.
- **Selection Type Classification**: Determines whether a gene is under **positive selection**, **negative selection**, or is **neutral** based on the dN/dS ratio.

## Installation
To run this project, you need to have Python installed along with the necessary dependencies.

1. Install **Biopython** and **CodonAlign** via pip:

```bash
pip install biopython
```

## Usage
### Step-by-Step Guide:
1. **Fetch GenBank data**: The program fetches genomic data for two viruses using their GenBank accession numbers (e.g., NC_045512.2 and PV009232.1).
2. **Gene Comparison**: The common genes between the two viral genomes are identified.
3. **Calculate dN/dS Ratio**: For each common gene, the program calculates the **dN/dS ratio** and classifies the gene into one of the following categories:
   - **Positive Selection**
   - **Negative Selection**
   - **Neutral Selection**
4. **Output**: The results are printed in a tabular format, showing the **gene name**, **dN**, **dS**, **dN/dS ratio**, and the type of selection acting on the gene.

### Example Usage:
To run the program, simply execute the script in Python:

```bash
python main.py
```

### Example Output:

```text
Synonymous positions:
{'ATA': 2, 'ATC': 2, 'ATT': 2, 'ATG': 0, 'ACA': 3, ...}

GenBank Data for NC_045512.2:
Sequence ID: NC_045512.2
Description: Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
Genome Length: 29903
Total number of genes: 23
Number of protein-coding genes: 12

GenBank Data for PV009232.1:
Sequence ID: PV009232.1
Description: Severe acute respiratory syndrome coronavirus 2 isolate SARS-CoV-2/human/USA/CA-LACPHL-AY09619/2024, complete genome
Genome Length: 29716
Total number of genes: 23
Number of protein-coding genes: 12

compare genes:
Common genes between NC_045512.2 and PV009232.1: {'ORF10', 'E', 'N', 'M', 'S', 'ORF3a', 'ORF1ab', 'ORF6', 'ORF7b', 'ORF8', 'ORF7a'}
Genes in NC_045512.2 but not in PV009232.1: set()
Genes in PV009232.1 but not in NC_045512.2: set()

Gene | dN | dS | dN/dS | Selection type
--------------------------------------------------
ORF10 | 0.056 | 0.018 | 3.111 | Positive selection
E | 0.045 | 0.025 | 1.800 | Positive selection
...
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing
Feel free to fork the repository, make changes, and submit pull requests. Contributions for improving functionality, bug fixes, and feature additions are always welcome.

## Contact
For any questions or suggestions, please contact the author at **shlomiasi1@gmail.com**.
