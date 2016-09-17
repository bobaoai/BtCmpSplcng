__author__ = 'tx'


def output_for_list(addr, aList, all_string=True):
    output = open(addr, "w")

    print "OUTPUT TEST"
    print aList[0]
    print len(aList)
    print "OUTPUT TEST END"
    if type(aList[0]) is list:
        # print "HAHA"
        for i in aList:

            # print "\t".join(i)
            if not all_string:
                output_list = list(i)
                for j in xrange(len(output_list)):
                    output_list[j] = str(output_list[j])
                i = output_list
            if i == []:
                output.write("\n")
            else:
                output.write("\t".join(i) + "\n")
    else:
        if type(aList[0]) == type("a"):
            for i in aList:
                output.write(i + "\n")
        else:
            try:
                for i in aList:
                    output.write(str(i) + "\n")
            except:
                print "Cannot print"


def output_for_fasta_table(address, aList):
    output = open(address, "w")
    print "OUTPUT FASTA TEST"
    print aList[0][0]
    print aList[0][1][0:10]
    print "OUTPUT FASTA TEST END"

    for name, seq in aList:
        # sequence = make_fasta_entry(name, seq)
        # print sequence
        output.write(make_fasta_entry(name, seq))

    print "There are " + str(len(aList)) + " lines outputed"
    output.close()


def output_for_fasta_dic(address, a_dic):
    output = open(address, "w")
    # print "OUTPUT FASTA TEST"
    #
    # print "OUTPUT FASTA TEST END"

    for name in a_dic.keys():
        a = make_fasta_entry(name, a_dic[name])
        # print a
        output.write(a)

    print "There are " + str(len(a_dic.keys())) + " lines outputed"
    output.close()


def make_fasta_entry(name, seq):
    return str(">" + name + "\n" + seq + "\n")
