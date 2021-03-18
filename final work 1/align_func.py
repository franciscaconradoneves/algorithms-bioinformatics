#common functions
def read_submat_file(filename):
    '''read substitution matrix from a file'''
    sm = {}
    f = open(filename, "r")
    line = f.readline()
    tokens = line.split("\t")
    ns = len(tokens)
    alphabet = []
    for i in range(0,ns):
        alphabet.append(tokens[i][0])
    for i in range(0,ns):
        line = f.readline()
        tokens = line.split("\t")
        for j in range(0,len(tokens)):
            k = alphabet[i]+alphabet[j]
            sm[k] = int(tokens[j])
    return sm  

def matrixSubs(alphabet, match, mismatch):
    ''' creates a substitution matrix given the values of match and mismatch'''
    dic = {}
    alphabet = alphabet.upper()
    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            if alphabet[i] == alphabet[j]:
                dic[alphabet[i]+alphabet[j]]= match
            else:
                dic[alphabet[i]+alphabet[j]]= mismatch
    return dic

def score_pos(c1,c2,sm, gap):
    '''score of a position(column)'''
    scoreValue = 0
    if c1 == "-" or c2 == "-":
        scoreValue += gap
    else: 
        scoreValue = sm[c1+c2]
    return scoreValue


def print_mat(mat):
    '''print in terminal row by row of matrix'''
    for i in range(0, len(mat)):
        print(mat[i])


def matrix(nrows, ncols):
    ''' creates a matrix nrows X ncols'''
    mt=[]
    for i in range(nrows):
        mt.append([0]*ncols)
    return mt

#global alignment

def max3t(v1,v2,v3):
    '''evaluate the integer or integers to fill in the T matrix'''
    maxT = max(v1,v2,v3)
    res = []
    if v1 == maxT:
        res.append(1)
    if v2 == maxT:
        res.append(2)
    if v3 == maxT:
        res.append(3)
    return res

def compare_pairwise_global_align(seqList, sm, g):
    '''return a matrix with all optimal scores found during the global alignment'''
    score_mat = [["-"]]
    for i in range (1,len(seqList)+1):
        score_mat[0].append(i)
        score_mat.append([i])
    #print(score_mat) 
    for firstSeq in range(0,len(seqList)):
        for secSeq in range(0, len(seqList)):
            #print("entrou no for ")
            res= global_align_multiple_solutions(seqList[firstSeq], seqList[secSeq], sm, g)
            S = res[0]
            score =S[len(seqList[firstSeq])][len(seqList[secSeq])]
            #print(score_mat[firstSeq])
            score_mat[firstSeq+1].append(score)
    return score_mat      
            

def global_align_multiple_solutions(seq1, seq2, sm, g):
    '''performs the global alignment between two sequences'''
    S = [[0]]
    T = [[0]]
    # initialize gaps in rows
    for j in range(1, len(seq2)+1):
        S[0].append(g * j)
        T[0].append(3)
        # initialize gaps in cols
    for i in range(1, len(seq1)+1):
        S.append([g * i])
        T.append([2])
        # apply the recurrence to fill the matrices
    for i in range(0, len(seq1)):    
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g)
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            S[i+1].append(max(s1, s2, s3))
            T[i+1].append(max3t(s1, s2, s3))
    return (S, T)


def recover_global_align_multiple_solutions(T,seq1, seq2):
    '''return all possible alignments between two sequences'''
    tmp_aligns = [["", "", len(seq1) , len(seq2) ]]
    final_aligns = []
    while tmp_aligns != []: 
        #print("temp_align", tmp_aligns)
        align = tmp_aligns.pop()
        i = align[2]
        j = align[3]
        if i == 0 and j == 0:
            final_aligns.append([align[0], align[1]])
            #print(align[0])
            #print(align[1]) 
        elif type(T[i][j]) == list:
            #print("entrou no elif")
            if i != 0 and j != 0:
                for t in T[i][j]:
                    #print("entrou no for")
                    if t == 1:
                        new_tmp_align = [seq1[i-1] + align[0], seq2[j-1] + align[1], i-1, j-1]
                    if t==3:
                        new_tmp_align = ["-" + align[0], seq2[j-1] + align[1], i, j-1]
                    if t==2:
                        new_tmp_align = [seq1[i-1] + align[0], "-" + align[1], i-1, j]        
                    tmp_aligns.append(new_tmp_align) 
                    #print(new_tmp_align)
        else:
            if T[i][j] == 1:
                new_tmp_align = [seq1[i-1] + align[0], seq2[j-1] + align[1], i-1, j-1]
            if T[i][j]==3:
                new_tmp_align = ["-" + align[0], seq2[j-1] + align[1], i, j-1]
            if T[i][j]==2:
                new_tmp_align = [seq1[i-1] + align[0], "-" + align[1], i-1, j]
            tmp_aligns.append(new_tmp_align)
            #print(new_tmp_align) 
    return final_aligns


#local alignment

def max3t_local(v1,v2,v3):
    '''evaluate the integer or integers to fill in the T matrix'''
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else: 
        if v2 >v3: return 2
        else: return 3

def local_align_multiple_solutions(seq1,seq2, sm, g):
    '''performs the local alignment between two sequences'''
    S= [[0]]
    T = [[0]]
    maxScore = 0
    for j in range(1,len(seq2)+1):
        S[0].append(0)
        T[0].append(0)
    for i in range(1,len(seq1)+1):
        S.append([0])
        T.append([0])
    for i in range(0,len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g)
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            b = max(s1,s2,s3)
            if b <= 0:
                S[i+1].append(0)
                T[i+1].append(0)
            else: 
                S[i+1].append(b)
                T[i+1].append(max3t_local(s1,s2,s3))
                if b > maxScore:
                    maxScore = b
    return (S,T,maxScore)

def compare_pairwise_local_align(seqList, sm, gap):
    '''return a matrix with all optimal scores found during the global alignment'''
    score_mat = [["-"]]
    for i in range(1, len(seqList)+1):
        score_mat[0].append(i)
        score_mat.append([i])
    print(score_mat)
    for seq1 in range(0, len(seqList)):
        for seq2 in range(0, len(seqList)):
            #print("entrou no for")
            res = local_align_multiple_solutions(seqList[seq1], seqList[seq2], sm, gap)
            score = res[2]
            #print(score)
            score_mat[seq1+1].append(score)
    return score_mat

def max_mat(mat):
    """finds the max cell in the matrix"""
    maxval = mat[0][0]
    maxrow = []
    maxcol = []
    for i in range(0,len(mat)):
        for j in range(0, len(mat[i])):
            if mat[i][j] >= maxval:
                maxval = mat[i][j]
    for i in range(0,len(mat)):
        for j in range(0, len(mat[i])):            
            if mat[i][j] == maxval:
                maxrow.append(i)
                maxcol.append(j)
    return (maxrow, maxcol)


def recover_local_align_multiple_solutions (S, T, seq1, seq2):
    """determine the cell with max score"""
    final_aligns = []
    m, n = max_mat(S)
    #print(m, n)
    for i in range(0,len(m)):
        tmp_aligns = [["","", m[i] ,n[i]]]
        #print(tmp_aligns)
        """terminates when finds a cell with zero"""
        while tmp_aligns != []:
            align = tmp_aligns.pop()
            i = align[2]
            j = align[3]
            if T[i][j] == 0:
                final_aligns.append([align[0], align[1]])
            elif type(T[i][j]) == list:
                if i!= 0 and j!= 0:
                #print("entrou no if é lista")
                    for t in T[i][j]:
                        #print("entrou no elif")
                        if T[i][j]==1:
                            new_tmp_align =[ seq1[i-1] + align[0], seq2[j-1] + align[1], i-1, j-1]
                        elif T[i][j] == 3:
                            new_tmp_align = ["-" + align[0], seq2[j-1] + align[1], i, j-1]
                        elif T[i][j] == 2:
                            new_tmp_align = [seq1[i-1] + align[0], "-" + align[1], i-1, j]
                        tmp_aligns.append(new_tmp_align)
                        #print(tmp_aligns)
            else:
                #print("entrou no else")
                if T[i][j]==1:
                    new_tmp_align =[ seq1[i-1] + align[0], seq2[j-1] + align[1], i-1, j-1]
                elif T[i][j] == 3:
                    new_tmp_align = ["-" + align[0], seq2[j-1] + align[1], i, j-1]
                elif T[i][j] == 2:
                    new_tmp_align = [seq1[i-1] + align[0], "-" + align[1], i-1, j]
                tmp_aligns.append(new_tmp_align)
                #print(tmp_aligns)
    return final_aligns




def test_align_global():
    sm = matrixSubs("ACTG", 1, -1)
    seq2= "GATTACA"
    seq1 = "GCATGCT"
    res = global_align_multiple_solutions(seq1, seq2, sm, -1)
    S = res[0]
    T = res[1]
    print("Score of optimal alignment: ", S[len(seq1)][len(seq2)])
    print_mat(S)
    print_mat(T)
    align =recover_global_align_multiple_solutions(T, seq1, seq2)
    print(align)



def test_align_local():
    sm = matrixSubs("ACTG", 1, -1)
    seq2= "GATTACA"
    seq1 = "GCATGCT"
    res =  local_align_multiple_solutions(seq1, seq2, sm, -1)
    S = res[0]
    T = res[1]
    print("Score of optimal alignment: ", res[2])
    print_mat(S)
    print()
    print_mat(T)
    align =recover_local_align_multiple_solutions(S,T, seq1, seq2)
    print(align)