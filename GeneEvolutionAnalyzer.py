import io
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from Bio.Align import PairwiseAligner
from Bio.Data import CodonTable

class GenBankHandler:
    def __init__(self, email):
        Entrez.email = email
    
    def fetch_data(self, accession_number):
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        data = handle.read()
        handle.close()
        handle = io.StringIO(data)
        return list(SeqIO.parse(handle, "genbank"))

class GenomeAnalysis:
    @staticmethod
    def get_genome_length(record):
        return len(record.seq)

    @staticmethod
    def count_genes(records):
        return sum(1 for record in records for feature in record.features if feature.type in {"gene", "CDS"})

    @staticmethod
    def count_protein_coding_genes(records):
        return sum(1 for record in records for feature in record.features if feature.type == "CDS" and "protein_id" in feature.qualifiers)
    
    @staticmethod
    def get_gene_names(records):
        return {feature.qualifiers.get("gene", ["unknown"])[0] for record in records for feature in record.features if feature.type in {"gene", "CDS"}}
    
    @staticmethod
    def compare_gene_sets(records_1, records_2, genome_1, genome_2):
        gene_set_1 = GenomeAnalysis.get_gene_names(records_1)
        gene_set_2 = GenomeAnalysis.get_gene_names(records_2)

        common_genes = gene_set_1.intersection(gene_set_2)
        unique_genes_1 = gene_set_1 - gene_set_2
        unique_genes_2 = gene_set_2 - gene_set_1
        
        print("\n===== Gene Comparison =====")
        print(f"Common genes between {genome_1} and {genome_2}: {common_genes}")
        print(f"Genes in {genome_1} but not in {genome_2}: {unique_genes_1}")
        print(f"Genes in {genome_2} but not in {genome_1}: {unique_genes_2}")
        
        return common_genes

class CodonAnalysis:
    @staticmethod
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

class SequenceProcessor:
    @staticmethod
    def remove_stop_codons(seq):
        stop_codons = {"TAA", "TAG", "TGA"}
        return "".join(seq[i:i+3] for i in range(0, len(seq) - 2, 3) if seq[i:i+3] not in stop_codons)
    
    @staticmethod
    def pad_sequence(seq):
        while len(seq) % 3 != 0:
            seq += "N"
        return seq

class SequenceAligner:
    @staticmethod
    def align_sequences(seq1, seq2):
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        alignments = list(aligner.align(seq1, seq2))
        best_alignment = alignments[0]
        return best_alignment.target[:len(best_alignment.query)], best_alignment.query[:len(best_alignment.target)]

class EvolutionaryAnalysis:
    def __init__(self):
        self.codon_table = CodonTable.unambiguous_dna_by_id[1]
    
    def calculate_dn_ds(self, common_genes, records_1, records_2):
        results = []
        for gene in common_genes:
            seq1, seq2 = self.extract_sequences(gene, records_1, records_2)
            if not seq1 or not seq2:
                continue
            seq1, seq2 = SequenceAligner.align_sequences(seq1, seq2)
            try:
                dN, dS = cal_dn_ds(CodonSeq(seq1), CodonSeq(seq2), codon_table=self.codon_table)
                dN_dS_ratio = dN / dS if dS > 0 else float('inf')
                selection_type = "Positive selection" if dN_dS_ratio > 1 else "Negative selection" if dN_dS_ratio < 1 else "Neutral selection"
                results.append({"Gene": gene, "dN": dN, "dS": dS, "dN/dS": dN_dS_ratio, "Selection type": selection_type})
            except Exception as e:
                print(f"Error processing gene {gene}: {e}")
        return results
    
    def extract_sequences(self, gene, records_1, records_2):
        seq1 = next((str(feature.location.extract(record).seq) for record in records_1 for feature in record.features if feature.type == "gene" and feature.qualifiers.get("gene", ["unknown"])[0] == gene), "")
        seq2 = next((str(feature.location.extract(record).seq) for record in records_2 for feature in record.features if feature.type == "gene" and feature.qualifiers.get("gene", ["unknown"])[0] == gene), "")
        seq1, seq2 = SequenceProcessor.remove_stop_codons(seq1), SequenceProcessor.remove_stop_codons(seq2)
        return SequenceProcessor.pad_sequence(seq1), SequenceProcessor.pad_sequence(seq2)

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

    synonymous_positions = {codon: CodonAnalysis.calculate_synonymous_positions(codon, gencode) for codon in gencode if gencode[codon] != '_'}
    print("===== Synonymous Positions =====")
    df = pd.DataFrame(list(synonymous_positions.items()), columns=["Codon", "Synonymous Positions"])
    print(df.to_string(index=False))

    handler = GenBankHandler("shlomiasi1@gmail.com")
    analysis = GenomeAnalysis()
    evolution = EvolutionaryAnalysis()
    
    accession_numbers = ["NC_045512.2", "PV009232.1"]
    records = {acc: handler.fetch_data(acc) for acc in accession_numbers}
    
    print("\n===== Gene Statistics =====")
    for acc, recs in records.items():
        print(f"GenBank Data for {acc}:")
        print(f"Total number of genes: {analysis.count_genes(recs)}")
        print(f"Number of protein-coding genes: {analysis.count_protein_coding_genes(recs)}")
    
    common_genes = analysis.compare_gene_sets(records[accession_numbers[0]], records[accession_numbers[1]], accession_numbers[0], accession_numbers[1])
    results = evolution.calculate_dn_ds(common_genes, records[accession_numbers[0]], records[accession_numbers[1]])
    
    print("\n===== Evolutionary Comparison =====")
    print("Gene       |  dN   |  dS   | dN/dS |  Selection type")
    print("-" * 50)
    for result in results:
        print(f"{result['Gene']:10} | {result['dN']:.3f} | {result['dS']:.3f} | {result['dN/dS']:.3f} | {result['Selection type']}")

if __name__ == "__main__":
    main()
