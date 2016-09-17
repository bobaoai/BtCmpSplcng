__author__ = 'tx'

# bb576

#1/30/2015 1st version for zea and musa alignment

#the psition is started with zero
def roastElement(fileaddress, columnlist, sep="\t", out=True):
    assert type(columnlist) == list and len(columnlist) >= 1
    file = open(fileaddress)

    returned_list = []
    line_acc = 0

    for i in file:
        line_acc += 1

        i = i[:-1]
        k = i.split(sep)

        saved_list = []
        for iter in columnlist:
            saved_list.append(k[iter])

        returned_list.append(saved_list)

    #check return
    if out:
        out_file = open(fileaddress + r".md", "w")
        for i in returned_list:
            out_file.write(sep.join(i) + "\n")
    return returned_list


if __name__ == "__main__":
    # address = r"/workdir/bb576/CML_500bp_S1_L001_R1_001.fastq"
    # out = r"/workdir/bb576/CML_500bp_S1_L001_R1_001.parse.fastq"
    # get_line()
    roastElement(r"D:\2015spring\Bioinfo\roast.chrom.10.msa.in.rates.full.elems", [1, 2])