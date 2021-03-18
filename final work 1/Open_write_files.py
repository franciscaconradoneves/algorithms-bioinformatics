def read_file(file_name):
    with open(file_name) as fh:
        #print("file open")
        lines = fh.readlines()
        seq = ""
        lst_Seq = []
        # print(lines)
        for line in range(1, len(lines)):
            #print("entrou no for ")
            if lines[line][:1] != ">":
                seq += lines[line].strip()
            else:
                lst_Seq.append(seq)
                seq = ""
                line += 1
        lst_Seq.append(seq)
    return lst_Seq


def write_file(seq, file_name):
    with open(file_name, "a") as fh:
        fh.write(seq + '\n')
    fh.close()