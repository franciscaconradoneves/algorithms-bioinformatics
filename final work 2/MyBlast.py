# -*- coding: utf-8 -*-

class MyBlast:
    '''
    Classe para matrizes de pontos
    '''

    def __init__(self, base_dados, w = 3):
        '''
        Construtor
        '''
        self.db = base_dados
        self.w = w
        self.map = None
    
    def rmSequenceDB(self, seq):
        """remove one sequence in DB"""
        res = []
        for i in self.db:
            if seq != i:
                if i not in res:   
                    res.append(i)
        self.db = res
    
    
    def build_map (self,query, w=3):
        res = {}
        for i in range(len(query[0][1])-w+1):
            subseq = query[0][1][i:i+w]
            if subseq in res:
                res[subseq].append(i)
            else:
                res[subseq] = [i]
        self.map = res
        #print(self.map)
        return res 

    
    def getHits (self, seq, query):
        res = [] # list of tuples
        for i in range(len(seq)-self.w+1):
            subseq = seq[i:i+self.w]
            if subseq in self.map:
                l = self.map[subseq]
                for ind in l:
                    res.append((ind,i))
        return res
        
    
    def extendsHit (self, seq, hit, query):
        stq, sts = hit[0], hit[1]
        ## move forward
        matfw = 0       
        k=0
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]: 
                matfw+=1
                bestk = k+1
            k += 1
        size = self.w + bestk
        ## move backwards
        k = 0
        matbw = 0   
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]: 
                matbw+=1
                bestk = k+1
            k+=1       
        size += bestk
        #print(stq-bestk, sts-bestk, size, self.w+matfw+matbw)
        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)
        
    def hitBestScore(self, seq, query):
        hits = self.getHits(seq, query)
        #print("Hits: ", hits)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        #print(best)
        return best
 
 
    def bestAlignment (self, query):
        self.build_map(query,self.w)
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            bestSeq = self.hitBestScore(self.db[k][1], query[0][1])
            if bestSeq != ():
                score = bestSeq[3] 
                if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res
