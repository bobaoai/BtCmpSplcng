from dataloading.build_file_for_fasta import read_seq_table
from dataloading.read_from_table import read_table
from dataloading.output_util import output_for_list

from dataloading.data_loading_util import head

__author__ = 'tx'


def blast_file_to_table(rawBlastnFile):
    """# input: raw blastn seq
    # work as filter: take the 100.00% identity seq out
    # print("Begin reading")
    # inputfile = open(rawBlastnFile)
    tag_count = 0
    # taglist = []
    tag_number = -1
    gene_name = []
    gene_value = 0
    tag_startposition = -1
    returnList = []
    # for i in inputfile:
    #     i = i[:-1]
    #     i = i.split("\t")
    #     identity = float(i[0])
    #
    #     tag_name = i[1]
    #     seq_name = i[2]
    #     start_position = i[3]
    #     gene_startpos = i[4]
    #     match_length = i[5]
    #     seq_length = i[6]
    #     #keywords includes:tag_number, matched_gene, start_position, match_length, gene_startpos
    #     tag_count = tag_count + 1
    # print("there are " + str(tag_count))
    # with open(IdentityFile, "w") as text_file:
    #     for i in tag_list:
    #         text_file.write(i.tag_number + " " + i.matched_gene+""+i.match_length+"\n")
    #     text_file.close()
    # print("end reading")"""
    aList = read_table(rawBlastnFile, sep="\t", head=False)
    return aList
    # return taglist


def filter_identity_ratio(aList, ratioThreshold=93.0):
    b_list = []
    for i in aList:
        # print float(i[0])
        if float(i[0]) < ratioThreshold and int(i[5])<400:
            continue
        b_list.append(i)
    print "Before filter we have alignments: ", len(aList)
    print "After filter we have alignments: ", len(b_list)
    return b_list



def count_gene_amount(aList, with_identity_ratio = True):
    # nameSet = [i[1] for i in aList]
    if with_identity_ratio:
        nameSet = [i[1] for i in aList]
    else:
        nameSet = [i[0] for i in aList]


    print "\n Name out of blast format List Check"
    for i in xrange(6):
        print nameSet[i]
    nameSet = list(set(nameSet))

    splitedNameSet = [i.split(".") for i in nameSet]

    orderList = [".".join(i[:2]) for i in splitedNameSet]

    # print "\n OrderList Check"
    # for i in xrange(6):
    #     print orderList[i]
    # print "\n splitedNameSet Check"
    # for i in xrange(6):
    #     print splitedNameSet[i]

    orderSet = [int(k[0]) for k in splitedNameSet]
    orderSet = list(set(orderSet))
    print "\nGenes that find at least one alignment", len(orderSet)

    # a1 = [str(i) + ".1" for i in orderSet]
    # a2 = [str(i) + ".2" for i in orderSet]
    # a1.extend(a2)

    misCount = 0
    misOrder = []

    for i in orderSet:
        if not str(i) + ".1" in orderList:
            misCount += 1
            misOrder.append(i)
            # print str(i) + ".1"
        if not str(i) + ".2" in orderList:
            misCount += 1
            misOrder.append(i)
            # print str(i) + ".2"

    # Since there are 3425 genes abstracted from the fasta
    # print "The result is for " + species
    print "Totally, " + str(misCount) + " gene has missed alignment\n", \
        "There are " + str(len(orderSet) - misCount) + " gene have both aligments", \
        "There are " + str(3425 - len(orderSet)) + " genes missing both alignments"
    removed_list = remove_missed_pair(aList, misOrder, misCount,with_identity_ratio=with_identity_ratio)
    print "\nAfter removed the alignments from genes with missed alignments.\nWe have alignments:  ", len(removed_list)
    return removed_list


def check_start_end(input_blastn_filtered_list):
    b_list = []
    print "file length", len(input_blastn_filtered_list)
    count = 0

    for i in input_blastn_filtered_list:
        if int(i[3]) > int(i[5]) or int(i[4]) > int(i[6]):
            # print i
            count += 1
            continue
        b_list.append(i)
    print "big",count
    return b_list

def remove_missed_pair(aList, misOrder, misCount, with_identity_ratio = True):
    b_list = list(aList)

    missedAmount = []

    for i in aList:
        # print i
        if with_identity_ratio:
            splited_name = i[1].split(".")
        else:
            splited_name = i[0].split(".")
        orderName = splited_name[0]
        # print orderName
        if int(orderName) in misOrder:
            missedAmount.append(orderName)
            b_list.remove(i)
            # print "\t".join(i)
    assert len(list(set(missedAmount))) == misCount, "MisCount not equal to missedAmount"
    return b_list
    # print max(g)
    # print len(x)
def remove_self_alignment(aList):
    bList = []
    for i in aList:
        qid_splited = i[1].split(".")
        if qid_splited[1]=="2":
            continue
        seqid_splited = i[2].split(".")
        if i[1]!= i[2]: bList.append(i)
    return bList

def remove_dif_order_alignment(aList):
    bList = []
    for i in aList:
        qid_splited = i[1].split(".")
        seqid_splited = i[2].split(".")
        if qid_splited[0] != seqid_splited[0]:
            # print i
            continue
        if i[1]!= i[2]: bList.append(i)
    return bList



def fasta_file_rebuilding_by_self_blast():
    global species
    aList = blast_file_to_table(rawBlastnFile=r"D:\2015spring\Bioinfo\VossiaCoelorachis\result_" + species + version+".txt")
    head(aList)

    print len(aList)
    self_alignments_removed = remove_self_alignment(aList)
    dif_alignment_removed = remove_dif_order_alignment(self_alignments_removed)
    ratio_filterd_self_alignment_removed = filter_identity_ratio(dif_alignment_removed, 0)
    final_gene_list = check_start_end(ratio_filterd_self_alignment_removed)

    gene_name_list = [[i[1],i[3],i[7]] for i in final_gene_list]
    print len(aList), len(self_alignments_removed), len(dif_alignment_removed),len(ratio_filterd_self_alignment_removed),len(gene_name_list)
    head(gene_name_list)
    geneFastaList = read_seq_table()
    output_for_list(r"D:\2015spring\Bioinfo\VossiaCoelorachis\gene_name_check",final_gene_list)
    outputAddr = r"D:\2015spring\Bioinfo\VossiaCoelorachis\HC_seq"+species+".txt"
    rebuild_fasta_by_gene_name_without_order(gene_name_list, geneFastaList, outputAddr)


def rebuild_fasta_by_gene_name_without_order(name_list_for_fasta, gene_cds_seq_dic, outputAddr):
    output = open(outputAddr, "w")
    # abstract_gene_name_fasta_title = []
    # both_exist = 0
    print "\nRebuild fasta"

    for i in name_list_for_fasta:

        qid_splited = i[0].split(".")
        name = ".".join(qid_splited[2:])
        sequence = gene_cds_seq_dic[name]
        start = int(i[1])
        length = int(i[2])
        end = start + length
        seq_required = sequence[start:end]
        print name, start, end, length
        name_string = "_".join([name,str(start),str(length)])
        output.write(">" + name_string+ "\n" + seq_required + "\n")




def blast_ana_pipeline():

    aList = blast_file_to_table(rawBlastnFile=r"D:\2015spring\Bioinfo\VossiaCoelorachis\result_" + species + version+".txt")

    bList = filter_identity_ratio(aList, 90)
    # filteredAlignmentFile = count_gene_amount(bList)
    # check_start_end(filteredAlignmentFile)

    output_for_list(addr="D:\\2015spring\\Bioinfo\\VossiaCoelorachis\\" + species + version+".filteredAlignments.txt",
               aList=bList)



def blast_ana_pipeline_for_coelorachis_6_22():

    aList = blast_file_to_table(rawBlastnFile=r"D:\2015spring\Bioinfo\Coelorachis\result_" + species + version+".txt")

    bList = filter_identity_ratio(aList, 80.0)
    filteredAlignmentFile = count_gene_amount(bList)
    # check_start_end(filteredAlignmentFile)

    output_for_list(addr="D:\\2015spring\\Bioinfo\\Coelorachis\\" + species + version+".filteredAlignments.txt",
               aList=bList)
if __name__ == "__main__":
    global species
    global version
    version = ".6.22"
    species = "Coelorachis.filtered"
    # species = "Sorghum"
    # species = "Vossia"
    # species = "self_p4"
    # fasta_file_rebuilding_by_self_blast()

    # for species in ["Vossia", "Coelorachis", "Sorghum", "CML247"]:
    #     blast_ana_pipeline()
    blast_ana_pipeline_for_coelorachis_6_22()
    # result_coelorachis.6.22.txt