
from MyAlign import MyAlign
from SubstMatrix import SubstMatrix

class PairwiseAlignment:

    def __init__(self, sm, g):
        self.g = g
        self.sm = sm
        self.S = None
        self.T = None
        self.seq1 = None    #list
        self.seq2 = None    #list
        
    def score_pos (self, c1, c2):
        if c1 == "-" or c2=="-":
            return self.g
        else:
            return self.sm[c1,c2]
        
    def score_alin (self, alin):
        res = 0
        
        for i in range(len(alin)):
            res += self.score_pos (alin[0][i], alin[1][i])
        return res
    
    def needleman_Wunsch (self, seq1, seq2):
        self.S = [[0]]
        self.T = [[0]]
        self.seq1 = seq1    #lista
        self.seq2 = seq2    #lista

        for j in range(1, len(seq2[1])+1):
            self.S[0].append(self.g * j)
            self.T[0].append(3)

        for i in range(1, len(seq1[1])+1):
            self.S.append([self.g * i])
            self.T.append([2])

        for i in range(0, len(seq1[1])):
            for j in range(len(seq2[1])):
                s1 = self.S[i][j] + self.score_pos (seq1[1][i], seq2[1][j])
                s2 = self.S[i][j+1] + self.g
                s3 = self.S[i+1][j] + self.g
                self.S[i+1].append(max(s1, s2, s3))
                self.T[i+1].append(max3t(s1, s2, s3))
        return self.S[len(seq1[1])][len(seq2[1])]
    

    def recover_align (self):
        res = []
        aux1 = [self.seq1[0],""]
        aux2 = [self.seq2[0],""]

        i = len(self.seq1[1])
        j = len(self.seq2[1])

        while i>0 or j>0:
            if self.T[i][j]==1:
                aux1[1] = self.seq1[1][i-1] + aux1[1]
                aux2[1] = self.seq2[1][j-1] + aux2[1]
                i -= 1
                j -= 1
            elif self.T[i][j] == 3:
                aux1[1] = "-" + aux1[1]
                aux2[1] = self.seq2[1][j-1] + aux2[1] 
                j -= 1
            else:
                aux1[1] = self.seq1[1][i-1] + aux1[1]
                aux2[1] = "-" + aux2[1]
                i -= 1
        res.append(aux1)
        res.append(aux2)
        print(res) 
        return MyAlign(res)

def max3t (v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3

def printMat (mat):
    for i in range(0, len(mat)):
        print(mat[i])


