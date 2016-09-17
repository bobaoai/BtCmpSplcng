from blast_result_analysis import blast_file_to_table, filter_identity_ratio
from dataloading.build_file_for_fasta import read_seq_table, checkRedundance, make_fasta_file_by_name_list
from dataloading.read_from_table import read_table, output_for_list

__author__ = 'tx'

def get_single_missed_gene(aList):
    nameSet = [i[1] for i in aList]

    print "\n Name out of blast format List Check"
    for i in xrange(6):
        print nameSet[i]
    nameSet = list(set(nameSet))

    splitedNameSet = [i.split(".") for i in nameSet]

    orderList = [".".join(i[:2]) for i in splitedNameSet]

    print "\n OrderList Check"
    for i in xrange(6):
        print orderList[i]
    print "\n splitedNameSet Check"
    for i in xrange(6):
        print splitedNameSet[i]

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
            misOrder.append(str(i) + ".1")
            # print str(i) + ".1"
        if not str(i) + ".2" in orderList:
            misCount += 1
            misOrder.append(str(i) + ".2")
            # print str(i) + ".2"

    # Since there are 3425 genes abstracted from the fasta
    # head(misOrder)
    print "missed gene's order" + str(misCount)
    return misOrder


def get_name_missed(misOrder, abstract_gene_list):
    global species
    # head(abstract_gene_list)
    split_gene_name_list = [i.split(".") for i in abstract_gene_list]
    # head(split_gene_name_list)
    order_and_name_list = [[".".join([i[0], i[1]]), i[2]] for i in split_gene_name_list]
    print "order_and_name_list"
    # head(order_and_name_list)
    single_missed_gene_name_list = []
    single_missed_gene_name_with_order_list = []
    for i in xrange(len(order_and_name_list)):
        if order_and_name_list[i][0] in misOrder:
            single_missed_gene_name_list.append(order_and_name_list[i][1])
            single_missed_gene_name_with_order_list.append(".".join(order_and_name_list[i]))
    # head(single_missed_gene_name_list)
    output_for_list("D:\\2015spring\\Bioinfo\\VossiaCoelorachis\\single_missed_gene." + species + ".txt",
               single_missed_gene_name_list)
    output_for_list("D:\\2015spring\\Bioinfo\\VossiaCoelorachis\\single_missed_gene_with_order." + species + ".txt",
               single_missed_gene_name_with_order_list)


def get_single_missed_gene_pipeline():
    global species


    geneAnnoList = read_table(r"D:\2015spring\Bioinfo\VossiaCoelorachis\sd01.txt", head=True)
    geneName = [x[:2] for x in geneAnnoList]
    # for i in geneName:
    # if i[0] == i[1]:
    # print "Yes"
    # print
    geneforFasta = checkRedundance(geneName)
    outputAddr = ""
    # geneFastaList = sort_fasta_to_seq_table(r"D:\2015spring\Bioinfo\VossiaCoelorachis\Zea_mays.AGPv3.23.cdna.all.fa")
    geneFastaList = read_seq_table()
    abstract_gene_with_order_list = make_fasta_file_by_name_list(geneforFasta, geneFastaList,
                                                                 r"D:\2015spring\Bioinfo\VossiaCoelorachis\mapFasta.txt")
    aList = blast_file_to_table(rawBlastnFile=r"D:\2015spring\Bioinfo\VossiaCoelorachis\result_" + species + ".6.8.txt")
    bList = filter_identity_ratio(aList)
    single_missed_gene_order_list = get_single_missed_gene(bList)
    get_name_missed(single_missed_gene_order_list, abstract_gene_with_order_list)