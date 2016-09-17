__author__ = 'tx'

def head(f, is_dic=False, title=""):
    if title != "":
        print title
    print "HEAD "
    k = 0
    if is_dic:
        print "Head for dic"
        line_ac = 0
        for i in f.keys():
            if f[i] is list:
                print i + ":"
                head(f[i])
            else:
                print i
                print f[i]
            line_ac += 1
            if line_ac == 10:
                return
    else:
        for i in f:
            print i
            if k > 10:
                break
            k += 1


def output_for_table(addr, aList, all_string = True):
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
                    output_list[j]= str(output_list[j])
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