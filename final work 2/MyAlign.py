class MyAlign:

    def __init__(self, lseqs, al_type = "protein"):
        self.listseqs = lseqs
        self.al_type = al_type
        self.testnumber = 0
    
    def __len__(self): # number of columns
        return len(self.listseqs[0])
    
    def __getitem__(self, n):
        if type(n) is tuple and len(n) ==2: 
            i, j = n
            return self.listseqs[i][j]
        elif type(n) is int: return self.listseqs[n]
        return None
    
    def __str__(self):
        res = " "
        for seq in self.listseqs:
            aux = " "
            aux += seq[0] + " " + seq[1]
            res += "\n" + aux 
        return res
    
    def num_seqs(self):
        return len(self.listseqs)
   
       
    def column (self, indice):
        res = []
        for k in range(len(self.listseqs)):
            res.append(self.listseqs[k][indice])
        return res
    
    def consensus(self):
        seq = []
        for i in self.listseqs:
            seq.append(i[1])
        cons = ["consensus",""]
        #print(seq[0])
        #print(seq[1])
        for i in range(491):    #491 tamanho das sequencias médias (com len nao estava a funcionar)
            cont = {}
            for k in range(len(seq)):
                c = seq[k][i]
                if c in cont:
                    cont[c] = cont[c] + 1
                else:
                    cont[c] = 1
            maximum = 0
            cmax = None
            for ke in cont.keys():
                if ke != "-" and cont[ke] > maximum:
                    maximum = cont[ke]
                    cmax = ke
            cons[1] = cons[1] + cmax
        #print(cons)
        return cons


    
    
        
    



##TEST##
'''
if __name__ == "__main__": 
    alig = MyAlign(["ATGA-A","AA-AT-"], "dna")
    print(alig)
    print(len(alig))
    print(alig.column(2))
    print(alig[1,1])
    print(alig[0])
    print(alig.consensus())

    alig2 = MyAlign(["VJKK","JRSK","VRSK"])
    print(alig2.richpolarbasic())
'''