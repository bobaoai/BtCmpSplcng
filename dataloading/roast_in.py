__author__ = 'tx'

# bb576
import os
# 1/30/2015 1st version for zea and musa alignment

#the psition is started with zero
# msa.in is fasta file
def roastIn(fileaddress, columnlist, out=True):
    print "Start roasr.in file loading"
    assert type(columnlist) == list and len(columnlist) >= 1
    number_genome = len(columnlist)
    file = open(fileaddress)

    returned_list = []
    line_acc = 0

    gen_seq = []
    gen_name = ""
    inlist = []
    #for storing seq
    for i in file:
        line_acc += 1
        i = i[:-1]
        # print len(i)
        if i.find(">") != -1:
            if gen_seq != "" and gen_name != "":
                if gen_name in columnlist:
                    returned_list.append([gen_name, gen_seq])
            gen_seq = ""
            gen_name = i[i.find(">") + 1:]
            #print gen_name
        elif gen_name not in columnlist:
            continue
        else:
            gen_seq += gen_seq + i
    if gen_name in columnlist:
        returned_list.append([gen_name, gen_seq])

    #check return
    if out:
        out_file = open(fileaddress + r".sel", "w")
        for i in returned_list:
            out_file.write(">" + i[0] + "\n")
            out_file.write(i[1] + "\n")
    # print "End roast.in file loading"
    return returned_list


def splitSpeciesSequencefromIn(fileaddress, outputname):
    print "Start roasr.in file loading"
    file = open(fileaddress)
    dirname = os.path.dirname(fileaddress)

    outputdir = os.path.join(dirname, outputname)
    if not os.path.exists(outputdir):
        os.mkdir(outputdir, 0754)

    returned_list = []

    gen_seq = []
    gen_name = ""

    line_acc = 0
    #for storing seq
    for i in file:
        line_acc += 1

        i = i[:-1]
        # print len(i)

        # new species name
        if i.find(">") != -1:

            #if list is not a new one
            if gen_seq != "" and gen_name != "":
                out_file = open(os.path.join(outputdir, gen_name), "w")
                out_file.write(gen_seq)
                # a new one
            gen_seq = ""
            gen_name = i[i.find(">") + 1:]
            # print gen_name

        else:
            gen_seq += gen_seq + i
    # the last one
    returned_list.append([gen_name, gen_seq])
    out_file = open(os.path.join(outputdir, gen_name), "w")
    out_file.write(gen_seq)
    # print "End roast.in file loading"


def loadSpeciesGenome(chromNum, speciesA, rootAddress=r"D:\2015spring\Bioinfo\\"):

    a = "\\" + "\\"
    address = a.join([rootAddress, "chrom." + str(chromNum), speciesA])
    file = open(address)
    return file.readline()


if __name__ == "__main__":
    # address = r"/workdir/bb576/CML_500bp_S1_L001_R1_001.fastq"
    # out = r"/workdir/bb576/CML_500bp_S1_L001_R1_001.parse.fastq"
    # # get_line()
    for i in range(1, 9):
        splitSpeciesSequencefromIn(r"D:\2015spring\Bioinfo\roast.chrom." + str(i) + ".msa.in", "chrom." + str(i))
        # a = []
        # for i in range(174, 266, 2):
        #     a.append(str(i))
        # print("".join(a))