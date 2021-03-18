
from PairwiseAlignment import PairwiseAlignment
from MyAlign import MyAlign
from SubstMatrix import SubstMatrix

class MultipleAlignment:

    def __init__(self, seqs, alignseq):
        self.seqs = seqs # list of MySeq objects
        self.alignpars = alignseq # PairwiseAlignment objects
    
    def num_seqs(self):
        return len(self.seqs)
    
    def add_seq_alignment (self, alignment, seq):
        res = []

        for i in range(len(alignment.listseqs)+1):
            res.append([])
            #print(res) criar as sublistas 

        # create consensus from give alignments
        #print(alignment)
        cons = alignment.consensus()
        self.alignpars.needleman_Wunsch(cons, seq)
        align2 = self.alignpars.recover_align()
        #print(type(align2))
        orig = 0
        
        for i in range(len(align2)):
            #print(align2[i,1])
            for j in range(len(align2[i,1])):
                #print(align2[i,1][j])
                if align2[i,1][j] == '-':
                    for k in range(len(alignment.listseqs)):
                        res[k] += "-"
                else:
                    for k in range(len(alignment.listseqs)):
                        res[k] += alignment[k,orig]
                    orig+=1
        #print(res)
        res[len(alignment.listseqs)] = align2.listseqs[1]
        return MyAlign(res)
    
    def align_consensus(self):
        self.alignpars.needleman_Wunsch(self.seqs[0], self.seqs[1])
        res = self.alignpars.recover_align()
        
        for i in range(2, len(self.seqs)):
            #print(self.seqs[i])
            res = self.add_seq_alignment(res, self.seqs[i])
        
        print("The Result from Multiple Alignment")
        print()
        print(res)
        return res
   

def printMat (mat):
    for i in range(0, len(mat)):
        print(mat[i])

