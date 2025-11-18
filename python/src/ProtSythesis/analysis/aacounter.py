class AACounter:

    aa_order = ("A","V","L","I","M","F","W","P","H","K","R","D","E","S","T","C","Y","N","Q","G","U")

    __slots__ = ("protein_id", "sequence", "seq_length",*aa_order)

    def __init__(self, protein_id, sequence):

        self.protein_id = protein_id
        self.sequence = sequence
        self.seq_length = len(sequence)

        self.A = 0
        self.V = 0
        self.L = 0
        self.I = 0
        self.M = 0
        self.F = 0
        self.W = 0
        self.P = 0 
        self.H = 0
        self.K = 0
        self.R = 0
        self.D = 0
        self.E = 0
        self.S = 0 
        self.T = 0
        self.C = 0
        self.Y = 0
        self.N = 0
        self.Q = 0
        self.G = 0
        self.U = 0
    
        self.count_aa()

    def count_aa(self):
        for aa in self.sequence:
            setattr(self, aa, getattr(self, aa) + 1)

    def get_normalized_counts(self):
        return [getattr(self,aa) / self.seq_length for aa in self.aa_order]
            



