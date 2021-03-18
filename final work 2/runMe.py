import re
import MyBlast
import MultipleAlignment
import MyAlign
import PairwiseAlignment
import SubstMatrix


def read_files(filename):
    headers = []
    seqs = []
    seq_aux = ""
    res = []
    #garantir que na primeira linha nada e adicionado ao seqs para nao ter conflitos depois
    flag = 0
    with open(filename) as fp:
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if flag != 0:
                    seqs.append(seq_aux)
                    seq_aux = ""
                header = re.search('\[([^]]+)\]', line)
                header = header.group()[1:-1]
                header = header.replace(" ", "_")   # substituicao do espaço por um underscore
                headers.append(header)
            else:
                flag = 1
                seq_aux += line
                #print(seq_aux)

        #garantir que a ultima seq_aux completa e adicionada na lista seqs 
        seqs.append(seq_aux)
        #print(seqs)

        for i in range(len(headers)):
            res.append([headers[i],seqs[i]])
    return res


def write_file(seq, file_name):
    with open(file_name, "a") as fh:
        fh.write(seq + '\n')
    fh.close()

def top10(base_dados,referencia):
    '''retorna a ref e as 10 sequencias mais proximas da ref '''
    res_blast = MyBlast.MyBlast(base_dados,11)
    res_blast.rmSequenceDB(referencia[0])
    res = referencia

    while len(res) < 11:
        a = res_blast.bestAlignment(referencia)
        seq_to_rm = res_blast.db[a[-1]]
        res.append(seq_to_rm)
        res_blast.rmSequenceDB(seq_to_rm)
        #print(res_blast.db)
    return res


def print_menu():       
    print (20 * "-" , "Welcome to Automatic Sequences Analysis Tool" , 20 * "-")
    print()
    print ("1. Menu Option to run Examples")
    print ("2. Menu Option to Start Analysis")
    print ("3. Exit")
    print (86 * "-")

loop=True      

while loop:          
    print_menu()   #Displays menu
    choice = int(input("Enter your choice [1-3]: "))
    if choice == 1:  
        print("Run Examples has been selected")
        print("To Blast analysis as reference we use the source.fasta file, as database we use the seqdump.txt file and")
        referencia = read_files("source.fasta")
        base_dados = read_files("seqdump.txt")
        top10Seq = top10(base_dados,referencia)     #resultado do blast top10 
        print("Results from Blast analysis were the first sequence is your reference ")
        print()
        print(top10Seq)
        print()
        print("To progessive multiple alignment we use Blosum62.mat as substitution matrix")
        #Alinhamento multiplo das 11 seqs
        #sm sera a default blosum62
        sm = SubstMatrix.SubstMatrix()
        sm.read_submat_file("blosum62.mat","\t")
        aseq = PairwiseAlignment.PairwiseAlignment(sm, -8)
        top10_align = MultipleAlignment.MultipleAlignment(top10Seq,aseq)
        ma_align = top10_align.align_consensus()

    elif choice == 2:
        print("Let's start the analysis ")
        referencia_File = input("Please choose a file with the reference: ")
        base_dados_File = input("Please choose a file with sequences to the database: ")
        referencia = read_files(referencia_File)
        base_dados = read_files(base_dados_File)
        top10Seq = top10(base_dados,referencia)     #resultado do blast top10 
        print("Results from Blast analysis were the first sequence is your reference ")
        print( )
        print(top10Seq)
        print( )
        print("Progressive Multiple Alignment")
        sm_File = input("Please choose a file with the substitution matrix: ")
        sm = SubstMatrix.SubstMatrix()
        sm.read_submat_file(sm_File,"\t")
        aseq = PairwiseAlignment.PairwiseAlignment(sm, -8)
        top10_align = MultipleAlignment.MultipleAlignment(top10Seq,aseq)
        ma_align = top10_align.align_consensus()

    elif choice == 3:
        print ("Exit has been selected")
        loop=False # This will make the while loop to end as not value of loop is set to False
    else:
        print(type(choice))
        input("Wrong option selection. Enter any key to try again..")








#verificacao nas especies duplicadas 
#a = ['Pongo abelii', 'MELSVLLFLALLTGLLLLLVQGHPNTHGRLPPGPRPLPLLGNLLQMDRRGLLKSFLRFREKYGDVFTVHLGPRPVVMLCGVEAIREALVDKAEAFSGRGKIAMVDPVFRGYGVIFANGNRWKVLRRFSVTTMRDFGMGKRSVEERIQEEAQCLIQELRKSKGALMDPTFLFQSITANIICSIVFGKRFHYQDQEFLKILNLFYQTFSLVSSVFGQLFELFSGFLKYFPGAHRQVYKNLQEINAYIGHSVEKHRETLDPSTPKDLIDTYLLHMEKEKSNAHSEFSHQNLTLNTLSLFFAGTETTSTTLRYGFLLMLKYPHVAERVYREIEQVIGPHRPPELHDRAKMPYTEAVIHEIQRFADLLPMGVPHIVTQHTSFRGYIIPKDTEVFLILSTALRDPHYFEKPDAFNPDHFLDASGALKKNEAFIPFSLGKRICLGEGIARAELFLFFTTILQNFSVASPVAPEDIDLTPQECGVGKIPPTYQIRFLPR']
#b = ['Pongo abelii', 'MELSVLLFLALLTGLLLLLVQRHPNTHGRLPPGPRPLPLLGNLLQMDRRGLLKSFLRFREKYGDVFTVHLGPRPVVVLCGVQAIREALVDKAEAFSGRGKIAIMDPVYQGYGVIFANGNRWKVLRRFSVTTMRDFGMGKQSVEERIQEEAQCLIEELQKSKGALMDPTFLFHSITANIICSIVFGKRFHYQDQEFLKMLNLFCQSFSLISSISSQLFELFSGFLKYFPGAHRQLYKNLQEINAYIGHSVEKHRETLDPSAPQDLIDTYLLHMEKEKSNPHSEFSHQNLIINTLSLFFAGTETTSTTLCYGFLLMLKYPHVAERVYKEIEQVVGPHCPPVLDDRAKMPYTEAVIHEIQRFADLLPMGVPHIVTQHTRFRGYIIPKDTEVFLILSTALRDPHYFEKPDAFNPDHFLDANGALKKNEAFIPFSLGKRICLGEGIARAELFLFFTTILQNFSVASPVAPEDIDLTPQECGVGKIPPTYQIRFLPH']
#c = ['Theropithecus gelada', 'MELSVLLFLALLTGLLLLLVQRHPNAHGRLPPGPRPLPLLGNLLQMDRRGLLRSFLRFREKYGDVFTVYLGPRPVVMLCGVEAIREALVDNAEAFSGRGKIAITDPVFQGYGVVFANGNRWKVLRRFSLTTMRDFGMGKRSVEERIQEEAQCLIEELRKSKGALVDPTFLFHSITANIICSIVFGKRFHYQDQEFLKILNLFYHTFSLASSMFGQLFELLSGFLKYFPGAHRQVYKNLQEINAYIGHSVEKHRETLDPSAPQDLIDSYLLQMEKEKSNPHSEFSHRNLIINTLSLFFAGTETTSTTLRYGFLLMLKYPHVAERIYKEIEQVIGPHRPPALDDRAKMPYTEAVIHEIQRFADLLPMGVPHIVTQQTSFRGYIIPKDTEVFPLLSTALHDPHYFEKPDTFNPDHFLDANGALKKNEAFIPFSLGRRICLGEGIARNELFLFFTTILQNFSVASPVAPEDIDLTPQESGVGKIPPKYQIRFLPR']
#d = ['Theropithecus gelada', 'MELSVLLFLALLTGLLLLLVQRHPNAHGRLPPGPRPLPLLGNLLQMDRRGLLRSFLRFREKYGDVFTVYLGPRPVVMLCGVEAIREALVDNAEAFSGRGKIAITDPVFQGYGVVFANGNRWKVLRRFSLTTMRDFGMGKRSVEERIQEEAQCLIEELRKSKGALVDPTFLFHSITANIICSIVFGKRFHYQDQEFLKILNLFYHTFSLASSMFGQLFELLSGFLKYFPGAHRQVYKNLQEINAYIGHSVEKHRETLDPSAPQDLIDSYLLQMEKEKSNPHSEFSHRNLIINTLSLFFAGTETTSTTLRYGFLLMLKYPHVAERIYKEIEQVIGPHRPPALDDRAKMPYTEAVIHEIQRFADLLPMGVPHIVTQQTSFRGYIIPKDTEVFPLLSTALHDPHYFEKPDTFNPDHFLDANGALKKNEAFIPFSLGRRICLGEGIARNELFLFFTTILQNFSVASPVAPEDIDLTPQESGVGKIPPTYQIRFLPR']
#a == b #False logo são genes diferentes 
#c == d #False logo são genes diferentes 