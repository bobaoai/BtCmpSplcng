from dataloading.read_from_table import read_table
from dataloading.data_loading_util import head, make_fasta_entry

__author__ = 'tx'




def checkRedundance(geneNameCluster):
    geneNameList = []
    existedGene = []
    order = 0
    for i in geneNameCluster:
        order += 1
        if "" in i: continue
        if "chr" in i[0] or "chr" in i[1]: continue
        new_gene = []
        for j in i:
            if j.find("||") != -1:
                j = j.split("||")
                new_gene.append(j)
            else:
                new_gene.append([j])
        existsIf = False
        for gene in new_gene[0]:
            if gene in existedGene:
                existsIf = True
        for gene in new_gene[1]:
            if gene in existedGene:
                existsIf = True
        if existsIf:
            existsIf = False
            continue
        kk = new_gene[0][:]
        for gene in new_gene[1]:
            kk.append(gene)
        # print kk
        alist = list(set(kk))
        # print alist, type(alist)
        for gene in alist:
            existedGene.append(gene)

        geneNameList.append(new_gene)
    print "geneNameList length: ", len(geneNameList)
    for i in xrange(6):
        print geneNameList[i]
    return geneNameList


def read_seq_table():
    addr = r"D:\2015spring\Bioinfo\VossiaCoelorachis\FastaTabFile.txt"
    input = open(addr)
    fastaList = {}
    for i in input:
        k = i[:-1].split("\t")
        fastaList[k[0]] = k[1]
    return fastaList


def make_fasta_file_by_name_list(name_list_for_fasta, gene_cds_seq_dic, outputAddr):
    output = open(outputAddr, "w")
    abstract_gene_name_fasta_title = []
    acc = 0
    both_exist = 0

    print "collect missed gene"
    for i in name_list_for_fasta:
        acc += 1
        exists1 = False
        exists2 = False
        dynamic_write_list = []

        for a1 in i[0]:
            if a1 in gene_cds_seq_dic.keys():
                exists1 = True
                dynamic_write_list.append(make_fasta_entry(str(acc) + ".1." + a1, gene_cds_seq_dic[a1]))
                # else:
                # print "Missed" + str(acc) + ".1." + a1
                # else:
                # print str(acc) + ".1." + a1
        for a2 in i[1]:
            if a2 in gene_cds_seq_dic.keys():
                exists2 = True
                dynamic_write_list.append(make_fasta_entry(str(acc) + ".2." + a2, gene_cds_seq_dic[a2]))
                # else:
                # print "Missed" + str(acc) + ".2." + a2
                # else:
                # print str(acc) + ".2." + a2
        if exists2 and exists1:
            both_exist += 1
            abstract_gene_name_fasta_title.extend(dynamic_write_list)
            for gene in dynamic_write_list:
                output.write(gene)

    print "end print missed gene"
    print "Gene waited for fasta in total: " + str(acc)
    print "Gene successfully abstrated for fasta in total: " + str(both_exist)
    successfully_abstract_gene_name_with_order_list = [i[1:i.find("\n")] for i in abstract_gene_name_fasta_title]
    return successfully_abstract_gene_name_with_order_list



def data_merge_pipeline():
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
    make_fasta_file_by_name_list(geneforFasta, geneFastaList, r"D:\2015spring\Bioinfo\VossiaCoelorachis\mapFasta.txt")


if __name__ == "__main__":
    data_merge_pipeline()