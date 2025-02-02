import io
from Bio import Entrez, SeqIO
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

def calculate_synonymous_positions(codon, gencode):
    synonymous_count = 0
    amino_acid = gencode[codon]  
    bases = ['A', 'C', 'G', 'T']   
    for position in range(3):
        for base in bases:
            mutated_codon = list(codon)
            mutated_codon[position] = base
            mutated_codon = ''.join(mutated_codon)
            if mutated_codon != codon and gencode.get(mutated_codon) == amino_acid:
                synonymous_count += 1
    return synonymous_count

def fetch_genbank_data(accession_number):
    Entrez.email = "shlomiasi1@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    data = handle.read()
    handle.close()
    handle = io.StringIO(data)
    records = list(SeqIO.parse(handle, "genbank"))  # Convert the iterator to a list
    return records

def Get_Genome_length(record):
    return len(record.seq)

def count_total_genes(records):
    total_genes = 0
    for record in records:
        for feature in record.features:
            if feature.type == "gene" or feature.type == "CDS":
                total_genes += 1
    return total_genes

def count_proteins(records):
    protein_coding_genes = 0
    for record in records:
        for feature in record.features:
            if feature.type == "CDS" and ("protein_id" in feature.qualifiers):
                protein_coding_genes += 1
    return protein_coding_genes

def get_gene_names(records):
    gene_names = set() 
    for record in records:
        for feature in record.features:
            if feature.type == "gene" or feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                gene_names.add(gene_name)
    return gene_names

def compare_gene_sets(gene_set_1, gene_set_2, genome_1, genome_2):
    common_genes = gene_set_1.intersection(gene_set_2)
    unique_genes_1 = gene_set_1 - gene_set_2
    unique_genes_2 = gene_set_2 - gene_set_1
    print("compare genes:")
    print(f"Common genes between {genome_1} and {genome_2}: {common_genes}")
    print(f"Genes in {genome_1} but not in {genome_2}: {unique_genes_1}")
    print(f"Genes in {genome_2} but not in {genome_1}: {unique_genes_2}")
    return common_genes

def calculate_dn_ds_for_common_genes(common_genes, records_1, records_2):
    results = []
    stop_codons = {"TAA", "TAG", "TGA"}  # קודוני עצירה
    for gene in common_genes:
        seq1 = ""
        seq2 = ""      
        for record in records_1:
            for feature in record.features:
                if feature.type == "gene" and feature.qualifiers.get("gene", ["unknown"])[0] == gene:
                    seq1 = feature.location.extract(record).seq
                    break
        for record in records_2:
            for feature in record.features:
                if feature.type == "gene" and feature.qualifiers.get("gene", ["unknown"])[0] == gene:
                    seq2 = feature.location.extract(record).seq
                    break

        # סינון קודוני עצירה
        seq1 = "".join([codon for codon in str(seq1).split() if codon not in stop_codons])
        seq2 = "".join([codon for codon in str(seq2).split() if codon not in stop_codons])
        
        # ודא שאורכי הרצפים תואמים
        if len(seq1) % 3 != 0 or len(seq2) % 3 != 0:
            print(f"Skipping {gene} because its sequence length is not divisible by 3.")
            continue
        
        if len(seq1) != len(seq2):
            print(f"Skipping {gene} because sequence lengths do not match.")
            continue
        
        codon_seq1 = CodonSeq(str(seq1))
        codon_seq2 = CodonSeq(str(seq2))     
        
        try:
            # חישוב dn ו-ds
            dN, dS = cal_dn_ds(codon_seq1, codon_seq2)
            dN_dS_ratio = float(dN / dS)
            
            # סוג הסלקציה
            if dN_dS_ratio > 1:
                selection_type = "Positive selection"
            elif dN_dS_ratio < 1:
                selection_type = "Negative selection"
            else:
                selection_type = "Neutral selection"
            
            results.append({
                "Gene": gene,
                "dN": dN,
                "dS": dS,
                "dN/dS": dN_dS_ratio,
                "Selection type": selection_type
            })
        except Exception as e:
            print(f"Error processing gene {gene}: {e}")
            continue

    return results


def main():
    gencode = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
    }
    synonymous_positions = {codon: calculate_synonymous_positions(codon, gencode) for codon in gencode if gencode[codon] != '_'}
    print("Synonymous positions:")
    print(synonymous_positions)
    accession_numbers = ["NC_045512.2", "PV009232.1"]
    gene_sets = {}
    for accession in accession_numbers:
        records = fetch_genbank_data(accession)  
        total_genes = count_total_genes(records)
        protein_coding_genes = count_proteins(records)
        gene_sets[accession] = get_gene_names(records)
        for record in records:
            print(f"\nGenBank Data for {accession}:")
            print(f"Sequence ID: {record.id}")
            print(f"Description: {record.description}")
            print(f"Genome Length: {Get_Genome_length(record)}")        
        print(f"Total number of genes: {total_genes}")
        print(f"Number of protein-coding genes: {protein_coding_genes}")
    common_genes = compare_gene_sets(gene_sets[accession_numbers[0]], gene_sets[accession_numbers[1]], accession_numbers[0], accession_numbers[1])
    results = calculate_dn_ds_for_common_genes(common_genes, fetch_genbank_data(accession_numbers[0]), fetch_genbank_data(accession_numbers[1]))
    print("Gene | dN | dS | dN/dS | Selection type")
    print("-" * 50)
    for result in results:
        print(f"{result['Gene']} | {result['dN']:.3f} | {result['dS']:.3f} | {result['dN/dS']:.3f} | {result['Selection type']}")

if __name__ == "__main__":
    main()