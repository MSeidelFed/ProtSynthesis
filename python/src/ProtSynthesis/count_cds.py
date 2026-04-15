import pyfastx
import numpy as np
import polars as pl 
 
codon_to_aa = {"TTT" : "F",
          "TTC" : "F",
          "TTA" : "L",
          "TTG" : "L",
          "TCT" : "S",
          "TCC" : "S",
          "TCA" : "S",
          "TCG" : "S",
          "TAT" : "Y",
          "TAC" : "Y",
          "TGT" : "C",
          "TGC" : "C",
          "TGG" : "W",
          "CTT" : "L",
          "CTC" : "L",
          "CTA" : "L",
          "CTG" : "L",
          "CCT" : "P",
          "CCC" : "P",
          "CCA" : "P",
          "CCG" : "P",
          "CAT" : "H",
          "CAC" : "H",
          "CAA" : "Q",
          "CAG" : "Q",
          "CGT" : "R",
          "CGC" : "R",
          "CGA" : "R",
          "CGG" : "R",
          "ATT" : "I",
          "ATC" : "I",
          "ATA" : "I",
          "ATG" : "M",
          "ACT" : "T",
          "ACC" : "T",
          "ACA" : "T",
          "ACG" : "T",
          "AAT" : "N",
          "AAC" : "N",
          "AAA" : "K",
          "AAG" : "K",
          "AGT" : "S",
          "AGC" : "S",
          "AGA" : "R",
          "AGG" : "R",
          "GTT" : "V",
          "GTC" : "V",
          "GTA" : "V",
          "GTG" : "V",
          "GCT" : "A",
          "GCC" : "A",
          "GCA" : "A",
          "GCG" : "A",
          "GAT" : "D",
          "GAC" : "D",
          "GAA" : "E",
          "GAG" : "E",
          "GGT" : "G",
          "GGC" : "G",
          "GGA" : "G",
          "GGG" : "G",
          "TGA" : "U"} 



class CDSCounter:

    """
    This class counts codons and amino acids and stores them in numpy arrays with a fix position
    The class can be converted into a row of a dataframe using the .to_df() method.

    gene_id: str (ensembl gene_id derived from fasta)
    transcript id: str (ensembl transcript_id derived from fasta)
    cds: str (cds derived from fasta)
    
    """

    codon_positions = {c:i for i,c in enumerate(codon_to_aa.keys())}
    aa_positions = {a:i for i,a in enumerate(set(codon_to_aa.values()))}


    def __init__(self, gene_id, transcript_id, cds):

        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.codon_counts = np.zeros(len(self.codon_positions.keys()))
        self.aa_counts = np.zeros(len(self.aa_positions.keys()))
    

        self.count_codons_and_aas(cds)

    def count_codons_and_aas(self, cds):    

        self.cds_length = len(cds)
        self.codon_count = 0
   
        for n in range(3,self.cds_length-3, 3):
            codon = cds[n:n+3]
            
            try:
                aa = codon_to_aa[codon]
                self.codon_counts[self.codon_positions[codon]] += 1
                self.aa_counts[self.aa_positions[aa]] += 1
                self.codon_count += 1
            except:
                print(f'For transcript {self.transcript_id} {codon} not found')


    def to_df(self):
        return [self.gene_id, self.transcript_id, self.codon_count] + list(self.codon_counts.tolist()) + list(self.aa_counts.tolist())

    @property
    def df_columns(self):
        return ['gene_id','transctipt_id', 'codon_count', *self.codon_positions.keys(), *self.aa_positions.keys()]


    @classmethod
    def from_seq(cls, seq):
        seq_info = {k:v for k,v in [x.split(':') for x in seq.description.split() if len(x.split(':')) == 2]}
        return cls(seq_info['gene'], seq.name, seq.seq)
    




def run(fasta_path):
    cds_counts = [CDSCounter.from_seq(s) for s in pyfastx.Fasta(fasta_path)]
    df_gene= pl.DataFrame([c.to_df() for c in cds_counts],schema = cds_counts[0].df_columns, orient="row")
    df_gene = df_gene.sort('codon_count').group_by('gene_id').last()

    return df_gene

