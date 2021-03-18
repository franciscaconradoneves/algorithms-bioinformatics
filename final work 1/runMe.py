import Open_write_files as OWF
import align_func as func

def print_menu():       ## Your menu design here
    print (40 * "-" , "MENU" , 40 * "-")
    print ("1. Menu Option to Perform Global Alignment")
    print ("2. Menu Option to Perform Local Alignment")
    print ("3. Menu Option to run Examples")
    print ("4. Exit")
    print (86 * "-")

loop=True      

while loop:          
    print_menu()   #Displays menu
    choice = int(input("Enter your choice [1-4]: "))
    if choice == 1:     
        print ("Global Alignment has been selected")
        SeqFile = input("Please choose a file with sequences to align: ")
        seqs = OWF.read_file(SeqFile)
        subs_mat_decision = input("Do you have a substitution matrix? [y/n] ")
        if subs_mat_decision == "y":
            Subs_mat = input("Enter you substitution matrix file: ")
            sm = func.read_submat_file(Subs_mat)
        else:
            alphabet = input("Enter the alphabet to use in the substitution matrix: ")
            match = int(input("Enter the value to match: "))
            mismatch = int(input("Enter the value to mismatch: "))
            sm = func.matrixSubs(alphabet,match, mismatch)
        gap = int(input("Please enter a value for gap: "))
        print()
        print("The list of sequences to align are: ")
        func.print_mat(seqs)
        print(20 * "*", "Possible solutions to global alignment", 20 * "*")
        for first in seqs:
            for second in seqs:
                print("SEQUENCE: ", first)
                print("SEQUENCE: ", second)
                res = func.global_align_multiple_solutions(first, second, sm, gap)
                S = res[0]
                T = res[1]
                print("Score matrix")
                func.print_mat(S)
                print()
                print("Trace-back matrix")
                func.print_mat(T)
                print()
                print("The score for optimal global alignment is: ", S[len(first)][len(second)])
                align = func.recover_global_align_multiple_solutions(T, first, second)    
                func.print_mat(align)
        table = func.compare_pairwise_global_align(seqs, sm, gap)
        print("global alignment optimal scores between all sequences in ", SeqFile)
        func.print_mat(table)
        OWF.write_file("Sequences aligned", "res.txt")
        for i in range(len(seqs)):
            OWF.write_file(seqs[i], "res.txt")
        OWF.write_file(" ", "res.txt")
        OWF.write_file("Global alignment optimal scores between all sequences", "res.txt")
        for i in table:
                OWF.write_file(str(i),"res.txt")
        OWF.write_file(" ", "res.txt")
        print("The result of this analysis has been saved in res.txt")

    elif choice==2:
        print ("Local alignment has been selected")
        SeqFile = input("Please choose a file with sequences to align: ")
        seqs = OWF.read_file(SeqFile)
        subs_mat_decision = input("Do you have a substitution matrix? [y/n] ")
        if subs_mat_decision == "y":
            Subs_mat = input("Enter you substitution matrix file: ")
            sm = func.read_submat_file(Subs_mat)
        else:
            alphabet = input("Enter the alphabet to use in the substitution matrix: ")
            match = int(input("Enter the value to match: "))
            mismatch = int(input("Enter the value to mismatch: "))
            sm = func.matrixSubs(alphabet,match, mismatch)
        gap = int(input("Please enter a value for gap: "))
        print()
        print("The list of sequences to align are: ")
        func.print_mat(seqs)
        print(20 * "*", "Possible solutions to Local alignment", 20 * "*")
        for first in seqs:
            for second in seqs:
                print("SEQUENCE: ", first)
                print("SEQUENCE: ", second)
                res = func.local_align_multiple_solutions(first, second, sm, gap)
                S = res[0]
                T = res[1]
                print("Score matrix")
                func.print_mat(S)
                print()
                print("Trace-back matrix")
                func.print_mat(T)
                print()
                print("The score for optimal local alignment is: ", res[2])
                align = func.recover_local_align_multiple_solutions(S, T, first, second)    
                func.print_mat(align)
                scoreRow, scoreCol = func.max_mat(S)
                if len(scoreRow) > 1:
                    print("There are multiple optimal alignments between ", first, second)
        table = func.compare_pairwise_local_align(seqs, sm, gap)
        print("local alignment optimal scores between all sequences in ", SeqFile)
        func.print_mat(table)
        OWF.write_file("Sequences aligned", "res.txt")
        for i in range(len(seqs)):
            OWF.write_file(seqs[i], "res.txt")
        OWF.write_file(" ", "res.txt")
        OWF.write_file("Local alignment optimal scores between all sequences", "res.txt")
        for i in table:
                OWF.write_file(str(i),"res.txt")
        OWF.write_file(" ", "res.txt")
        print("The result of this analysis has been saved in res.txt")

    elif choice==3:
        print ("Run Examples has been selected")
        print()
        print(25 * "*", " Example 1 ", 25 * "*")
        OWF.write_file("Example 1", "res.txt")
        seqs = OWF.read_file("protein_sequences.fas")
        gap = -3 
        sm = func.read_submat_file("blosum62.mat")
        print()
        print("Using the sequences present in the file protein_sequences.fas")
        print("Using a gap of -3 and Blosum62 substitution matrix")
        print()
        print("Possible solutions for Global Alignment")
        for firstSeq in seqs:
                for secSeq in seqs:
                    print("SEQUENCE: ", firstSeq)
                    print("SEQUENCE: ", secSeq)
                    res = func.global_align_multiple_solutions(firstSeq,secSeq, sm, gap)
                    S = res[0]
                    T = res[1]
                    #print("Score matrix")
                    #func.print_mat(S)
                    print()
                    #print("Trace-back matrix")
                    #func.print_mat(T)
                    print()
                    print("The score for the optimal global alignment is: ", S[len(firstSeq)][len(secSeq)])
                    align = func.recover_global_align_multiple_solutions(T, firstSeq, secSeq)
                    func.print_mat(align)
        table = func.compare_pairwise_global_align(seqs,sm, gap)
        func.print_mat(table)
        OWF.write_file("Sequences aligned", "res.txt")
        for i in range(len(seqs)):
            OWF.write_file(seqs[i], "res.txt")
        OWF.write_file(" ", "res.txt")
        OWF.write_file("Global alignment optimal scores between all sequences", "res.txt")
        for i in table:
                OWF.write_file(str(i),"res.txt")
        OWF.write_file(" ", "res.txt")
        print("The result of this analysis has been saved in res.txt")
        print()
        print("Possible solutions for Local alignment")
        for firstSeq in seqs:
            for secSeq in seqs:
                print("SEQUENCE: ", firstSeq)
                print("SEQUENCE: ", secSeq)
                res = func.local_align_multiple_solutions(firstSeq, secSeq, sm, -8)
                print("Score of optimal alignment:", res[2])
                S = res[0]
                T = res[1]
                align = func.recover_local_align_multiple_solutions(S,T,firstSeq,secSeq)
                print(align)
                print()
                scoreRow, scoreCol = func.max_mat(S)
                if len(scoreRow) > 1:
                    print("There are multiple optimal alignments between ", firstSeq, secSeq)
        table = func.compare_pairwise_local_align(seqs, sm,gap)
        func.print_mat(table)
        OWF.write_file("Sequences aligned", "res.txt")
        for i in range(len(seqs)):
            OWF.write_file(seqs[i], "res.txt")
        OWF.write_file(" ", "res.txt")
        OWF.write_file("Local alignment optimal scores between all sequences", "res.txt")
        for i in table:
                OWF.write_file(str(i),"res.txt")
        print("The result of this analysis has been saved in res.txt")
        print()
        print()
        print(25 * "*", " Example 2 ", 25 * "*")
        OWF.write_file("Example 2", "res.txt")
        seq1 = "GATTACA"
        seq2 = "GCATGCT"
        gap = -1
        sm = func.matrixSubs("ATCG", 1, -1)
        print()
        print("Using the following sequences: 'GATTACA' and 'GCATGCT'")
        print("Using a gap of -1 and substitution matrix created using the values for match 1 and mismatch -1")
        print()
        print("Possible solutions for Global Alignment")
        res = func.global_align_multiple_solutions(seq1,seq2, sm, gap)
        S = res[0]
        T = res[1]
        print("Score matrix")
        func.print_mat(S)
        print()
        print("Trace-back matrix")
        func.print_mat(T)
        print()
        print("The score for the optimal global alignment is: ", S[len(seq1)][len(seq2)])
        align = func.recover_global_align_multiple_solutions(T, seq1, seq2)
        func.print_mat(align)
        OWF.write_file("Sequences aligned", "res.txt")
        OWF.write_file(seq1, "res.txt")
        OWF.write_file(seq2, "res.txt")
        OWF.write_file(" ", "res.txt")
        print("The result of this analysis has been saved in res.txt")
        OWF.write_file("Global alignment optimal scores between all sequences", "res.txt")
        OWF.write_file(str(S[len(seq1)][len(seq2)]),"res.txt")
        OWF.write_file(" ", "res.txt")
        print("The result of this analysis has been saved in res.txt")
        print()
        print("Possible solutions for Local alignment")
        res = func.local_align_multiple_solutions(seq1, seq2, sm, -8)
        print("Score of optimal alignment:", res[2])
        S = res[0]
        T = res[1]
        align = func.recover_local_align_multiple_solutions(S,T,seq1,seq2)
        print(align)
        print()
        OWF.write_file("Sequences aligned", "res.txt")
        OWF.write_file(seq1, "res.txt")
        OWF.write_file(seq2, "res.txt")
        OWF.write_file(" ", "res.txt")
        OWF.write_file("Local alignment optimal score between the two sequences", "res.txt")
        OWF.write_file(str(res[2]),"res.txt")
        print("The result of this analysis has been saved in res.txt")

    elif choice==4:
        print ("Exit has been selected")
        loop=False # This will make the while loop to end as not value of loop is set to False
    else:
        print(type(choice))
        input("Wrong option selection. Enter any key to try again..")


