from dataloading.data_acquisite_multialignment import get_seq_from_sorghum_helper
from dataloading.read_from_fasta import *
from dataloading.read_from_table import *
from dataloading.remove_redundant_scaffold_from_incomplete_genome import GeneVote, SegmentWithScaffoldInformation
# build dic for 2 column
from dataloading.output_util import *
from os import path

__author__ = 'tx'


# file = open(r"D:\2015spring\Bioinfo\Coelorachis\Source data\coelorachis.lines.fasta")
# outfileaddr = r"D:\2015spring\Bioinfo\Coelorachis\Source data\coelorachis.lines.fasta"
#
# for orderi, i in enumerate(file):
#     if orderi % 250000 == 0:
#         output_file = open(outfileaddr+str(orderi),"w")
#     output_file.write(i)


def blast_as_ana_pipe():
    # load_pep_file()
    # load_coelorachis()
    load_sorghum()



def load_sorghum():
    sorghum_file = read_fasta_as_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\Sorghum_3063pair.fasta")
    # print len(coelorachis_file)
    merge_helper_sorghum(sorghum_file)


from Bio import SeqIO


def merge_helper_sorghum(the_list):
    # for index, record in enumerate(SeqIO.parse(open("ls_orchid.gbk"), "genbank")):
    # print "index %i, ID = %s, length %i, with %i features" \
    #       % (index, record.id, len(record.seq), len(record.features))
    genome_addr = r"Sorghum_bicolor.Sorbi1.27.dna.toplevel.fa"
    a_list = []
    handle = open(path.join(root_address, genome_addr), "rU")
    for record in SeqIO.parse(handle, "fasta"):
        print record.id
        name_string = str(record.id)
        name_splited = name_string.split(" dna")
        name_splited = name_splited[0][:]
        print name_splited
        if str(name_splited) not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]:
            continue
        else:
            a_list.append([name_splited, str(record.seq)])
    for i in a_list:
        print i[0]
    print "end"

    # a_list = read_fasta_as_table(path.join(root_address, genome_addr))
    # print
    # get_seq_from_sorghum_helper(a_list)
    # for i in a_list:
    #     if i[0]  in ["1","2","3","4","5","6","7","8","9","10"]:
    #         print i[0]
    global genome_dic
    genome_dic = {}
    # multiseq_dic = {}
    for i in a_list:
        genome_dic[i[0]] = i[1]
    # seperate to index, order
    # by using dic
    paird_list = []
    for i in the_list:
        title_part = seperate_name(i[0])
        # print index
        title_part.append(i[1])
        paird_list.append(title_part)
    # head(paird_list)
    gene_dic = {}
    for i in paird_list:
        if i[0] not in gene_dic.keys():
            gene_dic[i[0]] = [i]
        else:
            gene_dic[i[0]].append(i)
    seq_return_list = []
    for i in gene_dic.keys():
        temp_class = sorghum_pair(gene_dic[i])
        temp_class.get_final_seq()
        if temp_class.get_final_seq():
            seq_return_list.append(temp_class.get_final_seq())
    seq_return_table = []
    for i in seq_return_list:
        # print i[0:6]
        name = ":".join(i[0:6])
        # print name
        seq_return_table.append([int(i[0]), [name, i[6]]])
    seq_return_table = sorted(seq_return_table, key=itemgetter(0))
    seq_list_sorted = [i[1] for i in seq_return_table]
    head(seq_list_sorted)
    print len(seq_return_list)
    output_for_fasta_table(path.join(root_address, "sorghum_3063pair_1110.txt"), seq_list_sorted)
    # output_for_fasta_table()
    return seq_return_table


def load_coelorachis():
    coelorachis_file = read_fasta_as_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\Coelorachis_3063pair.fasta")
    # print len(coelorachis_file)
    merge_helper(coelorachis_file)


def merge_helper(the_list):
    # genome_addr = r"Coelorachis.a.lines.fasta"
    # a_list = read_fasta_as_table(path.join(root_address, genome_addr))
    # get_seq_from_sorghum_helper(a_list)
    # global genome_dic
    # genome_dic = {}
    # multiseq_dic = {}
    # for i in a_list:
    #     genome_dic[i[0]] = i[1]
    # seperate to index, order
    # by using dic
    paird_list = []
    for i in the_list:
        title_part = seperate_name(i[0])
        # print index
        title_part.append(i[1])
        paird_list.append(title_part)
    # head(paird_list)
    gene_dic = {}
    for i in paird_list:
        if i[0] not in gene_dic.keys():
            gene_dic[i[0]] = [i]
        else:
            gene_dic[i[0]].append(i)
    seq_return_list = []
    for i in gene_dic.keys():
        temp_class = gene_pair(gene_dic[i])
        seq_return_list.append(temp_class.get_final_seq())
    seq_return_table = []
    for i in seq_return_list:
        # print i[0:6]
        name = ":".join(i[0:6])
        # print name
        seq_return_table.append([int(i[0]), [name, i[6]]])
    seq_return_table = sorted(seq_return_table, key=itemgetter(0))
    seq_list_sorted = [i[1] for i in seq_return_table]
    head(seq_list_sorted)
    print len(seq_return_list)
    output_for_fasta_table(path.join(root_address, "coelorachis_3063pair_1110.txt"), seq_list_sorted)
    # output_for_fasta_table()
    return seq_return_table


class gene_pair():
    def __init__(self, alist):
        self.basic_list = alist
        # self.get_final_seq()
        #

    # 7', '2', 'GRMZM2G082976', 'flattened_line_16396', '6200', '10901'
    def get_final_seq(self):
        # name_set = set([i[3] for i in self.basic_list])
        scaffold_name_dic = {}
        for i in self.basic_list:

            if i[3] not in scaffold_name_dic.keys():
                scaffold_name_dic[i[3]] = [i]
            else:
                scaffold_name_dic[i[3]].append(i)

        for i in scaffold_name_dic.keys():
            for orderj, j in enumerate(scaffold_name_dic[i]):
                scaffold_name_dic[i][orderj].append(len(scaffold_name_dic[i]))
                scaffold_name_dic[i][orderj][4] = int(j[4])
                scaffold_name_dic[i][orderj][5] = int(j[5])
                print scaffold_name_dic[i][orderj][5] - scaffold_name_dic[i][orderj][4] + 1, len(
                    scaffold_name_dic[i][orderj][6])
                assert scaffold_name_dic[i][orderj][5] - scaffold_name_dic[i][orderj][4] + 1 == len(
                    scaffold_name_dic[i][orderj][6]), ("start_end wrong", scaffold_name_dic[i][orderj][:6])

            if len(scaffold_name_dic[i]) == 1:
                print scaffold_name_dic[i][0][:6]
            else:
                # 7', '2', 'GRMZM2G082976', 'flattened_line_16396', '6200', '10901', start at 6200 and ends at 10902 not 10901

                sorted_list = sorted(scaffold_name_dic[i], key=itemgetter(4))
                for j in sorted_list:
                    print "sorted"
                    print j[4]
                final_seq = ""
                iteri = 0
                while len(sorted_list) > 1:
                    if sorted_list[iteri][5] < sorted_list[iteri + 1][4]:
                        # iteri += 1
                        print "check!"

                        print "end check"
                        print sorted_list[iteri][:6]
                        print sorted_list[iteri + 1][:6]
                        astring1 = sorted_list[iteri][6]
                        lst_end1 = sorted_list[iteri][5]
                        lst_start2 = sorted_list[iteri + 1][4]
                        print lst_start2 - lst_end1
                        gap_seq = "N" * (lst_start2 - lst_end1 - 1)

                        astring2 = sorted_list[iteri + 1][6][:]
                        # print len(astring2), gap
                        # assert gap == len(astring2), "scissor mistake"
                        sorted_list[iteri][6] = astring1 + gap_seq + astring2
                        sorted_list[iteri][5] = sorted_list[iteri + 1][5]

                        print len(astring1), len(gap_seq), len(astring2)
                        print len(sorted_list[iteri][6]), sorted_list[iteri][5] - sorted_list[iteri][
                            4] + 1, "cancatenated with gap seq"
                        sorted_list.pop(iteri + 1)
                    else:
                        if sorted_list[iteri][5] >= sorted_list[iteri + 1][5]:
                            sorted_list.pop(iteri + 1)
                        else:
                            print sorted_list[iteri][:6]
                            print sorted_list[iteri + 1][:6]
                            astring1 = sorted_list[iteri][6]
                            lst_end1 = sorted_list[iteri][5]
                            lst_end2 = sorted_list[iteri + 1][5]
                            gap = lst_end2 - lst_end1

                            astring2 = sorted_list[iteri + 1][6][-gap:]
                            print len(astring2), gap
                            assert gap == len(astring2), "scissor mistake"

                            sorted_list[iteri][6] = astring1 + astring2
                            sorted_list[iteri][5] = sorted_list[iteri + 1][5]
                            sorted_list.pop(iteri + 1)
                            print len(sorted_list[iteri][6]), sorted_list[iteri][5] - sorted_list[iteri][
                                4] + 1, "cancatenated"
                scaffold_name_dic[i] = sorted_list

                # iteri += 1
        comp_list = []
        for i in scaffold_name_dic.keys():
            start = scaffold_name_dic[i][0][4]
            end = scaffold_name_dic[i][0][5]
            scaffold_name_dic[i][0].append(end - start + 1)
            assert len(scaffold_name_dic[i]) == 1
            comp_list.append(scaffold_name_dic[i][0])
        # head(comp_list)
        # fori in8
        comp_list = sorted(comp_list, key=itemgetter(8), reverse=False)
        for i in comp_list:
            print i[:6], i[7], i[8]
        print "finish " + comp_list[0][0]
        assert comp_list[0][5] - comp_list[0][4] + 1 == len(comp_list[0][6])

        comp_list[0][4] = str(comp_list[0][4])
        comp_list[0][5] = str(comp_list[0][5])
        comp_list[0][2] = "coelorachis"
        return comp_list[0]



        # return []


class sorghum_pair(gene_pair):
    def __init__(self, alist):
        self.basic_list = alist
        # self.get_final_seq()
        #

    # 7', '2', 'GRMZM2G082976', 'flattened_line_16396', '6200', '10901'
    def get_final_seq(self):
        # name_set = set([i[3] for i in self.basic_list])
        scaffold_name_dic = {}
        for i in self.basic_list:

            if i[3] not in scaffold_name_dic.keys():
                scaffold_name_dic[i[3]] = [i]
            else:
                scaffold_name_dic[i[3]].append(i)

        for i in scaffold_name_dic.keys():
            for orderj, j in enumerate(scaffold_name_dic[i]):
                scaffold_name_dic[i][orderj].append(len(scaffold_name_dic[i]))
                scaffold_name_dic[i][orderj][4] = int(j[4])
                scaffold_name_dic[i][orderj][5] = int(j[5])
                start = scaffold_name_dic[i][orderj][4]
                end = scaffold_name_dic[i][orderj][5] + 1
                if scaffold_name_dic[i][orderj][5] - scaffold_name_dic[i][orderj][4] + 1 == len(
                        scaffold_name_dic[i][orderj][6]):

                    anewstring = genome_dic[scaffold_name_dic[i][orderj][3]][start:end]

                    if scaffold_name_dic[i][orderj][6] != anewstring:
                        print "outch!!"
                        print scaffold_name_dic[i][orderj][:6]
                        scaffold_name_dic[i][orderj][6] = anewstring
                else:
                    print "woo"
                    if len(genome_dic[scaffold_name_dic[i][orderj][3]]) < scaffold_name_dic[i][orderj][5]:
                        scaffold_name_dic[i][orderj][5] = len(genome_dic[scaffold_name_dic[i][orderj][3]]) - 1
                    print len(genome_dic[scaffold_name_dic[i][orderj][3]])
                    scaffold_name_dic[i][orderj][6] = genome_dic[scaffold_name_dic[i][orderj][3]][start:]

                assert scaffold_name_dic[i][orderj][5] - scaffold_name_dic[i][orderj][4] + 1 == len(
                    scaffold_name_dic[i][orderj][6]), ("start_end wrong", scaffold_name_dic[i][orderj][:6])

            if len(scaffold_name_dic[i]) == 1:
                print scaffold_name_dic[i][0][:6]
            else:
                # 7', '2', 'GRMZM2G082976', 'flattened_line_16396', '6200', '10901', start at 6200 and ends at 10902 not 10901

                sorted_list = sorted(scaffold_name_dic[i], key=itemgetter(4))
                for j in sorted_list:
                    print "sorted"
                    print j[4]
                final_seq = ""
                iteri = 0
                while len(sorted_list) > 1:
                    if sorted_list[iteri][5] < sorted_list[iteri + 1][4]:
                        # iteri += 1
                        print "check!"

                        print "end check"
                        if sorted_list[iteri][5] < sorted_list[iteri + 1][4]-100000:
                            return []
                            print "huge gap"
                        print sorted_list[iteri][:6]
                        print sorted_list[iteri + 1][:6]
                        astring1 = sorted_list[iteri][6]
                        lst_end1 = sorted_list[iteri][5]
                        lst_start2 = sorted_list[iteri + 1][4]
                        print lst_start2 - lst_end1
                        gap_seq = "N" * (lst_start2 - lst_end1 - 1)

                        astring2 = sorted_list[iteri + 1][6][:]
                        # print len(astring2), gap
                        # assert gap == len(astring2), "scissor mistake"
                        sorted_list[iteri][6] = astring1 + gap_seq + astring2
                        sorted_list[iteri][5] = sorted_list[iteri + 1][5]

                        print len(astring1), len(gap_seq), len(astring2)
                        print len(sorted_list[iteri][6]), sorted_list[iteri][5] - sorted_list[iteri][
                            4] + 1, "cancatenated with gap seq"
                        sorted_list.pop(iteri + 1)
                    else:
                        if sorted_list[iteri][5] >= sorted_list[iteri + 1][5]:
                            sorted_list.pop(iteri + 1)
                        else:
                            print sorted_list[iteri][:6]
                            print sorted_list[iteri + 1][:6]
                            astring1 = sorted_list[iteri][6]
                            lst_end1 = sorted_list[iteri][5]
                            lst_end2 = sorted_list[iteri + 1][5]
                            gap = lst_end2 - lst_end1

                            astring2 = sorted_list[iteri + 1][6][-gap:]
                            print len(astring2), gap
                            assert gap == len(astring2), "scissor mistake"

                            sorted_list[iteri][6] = astring1 + astring2
                            sorted_list[iteri][5] = sorted_list[iteri + 1][5]
                            sorted_list.pop(iteri + 1)
                            print len(sorted_list[iteri][6]), sorted_list[iteri][5] - sorted_list[iteri][
                                4] + 1, "cancatenated"
                scaffold_name_dic[i] = sorted_list

        # # iteri += 1
        comp_list = []
        for i in scaffold_name_dic.keys():
            start = scaffold_name_dic[i][0][4]
            end = scaffold_name_dic[i][0][5]
            scaffold_name_dic[i][0].append(end - start + 1)
            assert len(scaffold_name_dic[i]) == 1
            comp_list.append(scaffold_name_dic[i][0])
        # head(comp_list)
        # fori in8
        comp_list = sorted(comp_list, key=itemgetter(8), reverse=False)
        for i in comp_list:
            print i[:6], i[7], i[8]
        print "finish " + comp_list[0][0]
        # assert comp_list[0][5] - comp_list[0][4] + 1 == len(comp_list[0][6])

        comp_list[0][4] = str(comp_list[0][4])
        comp_list[0][5] = str(comp_list[0][5])
        comp_list[0][2] = "sorghum"
        # if len(comp_list[0][6])>100000:
        return comp_list[0]



        # return []
def sorghum_coelorachis():
    sorghum_file = read_fasta_as_table(path.join(root_address, "sorghum_3063pair_1110.txt"))
    coelorachis_file = read_fasta_as_table(path.join(root_address, "coelorachis_3063pair_1110.txt"))
    fasta_file = r"D:\2015spring\Bioinfo\Coelorachis\Source data\Zea_mays.AGPv3.27.3325genes_seq.fasta"
    fasta_table = read_fasta_as_table(fasta_file)
    fasta_dic = {}
    fasta_list = []
    for i in fasta_table:
        # print i[0:6]
        name = i[0].split(".")[0]
        order = i[0].split(".")[1]
        a_name = ":".join(i[0].split("."))
        # print name
        # if name not in fasta_dic.keys():
        #     fasta_dic[name] = {"1":"", "2":""}
        #     fasta_dic[name][order] = [a_name, i[1]]
        # elif fasta_dic[name][order] == "":
        #     fasta_dic[name][order] = [a_name, i[1]]
        # # elif len(fasta_dic[name][order][1]) >
        #


        fasta_list.append([int(name), [a_name, i[1]]])
    # fasta_list = sorted(fasta_list, key=itemgetter(0))
    sorghum_list = []
    for i in sorghum_file:
        # print i[0:6]
        name = i[0].split(":")[0]
        # print name
        sorghum_list.append([int(name), i])
    sorghum_list = sorted(sorghum_list, key=itemgetter(0))
    coelorachis_list = []
    for i in coelorachis_file:
        # print i[0:6]
        name = i[0].split(":")[0]
        # print name
        coelorachis_list.append([int(name), i])
    coelorachis_list = sorted(coelorachis_list, key=itemgetter(0))
    coelorachis_set = set([i[0]for i in coelorachis_list])
    head(coelorachis_set)
    sorghum_set = set([i[0]for i in sorghum_list])
    head(sorghum_set)
    rebuilt_fasta_list = []
    new_coelorachis_list = []
    new_coelorachis_set = []
    for i in coelorachis_list:
        if i[0] in sorghum_set:
            new_coelorachis_set.append(i[0])
            new_coelorachis_list.append(i[1])
    print len(new_coelorachis_set)
    for i in fasta_list:
        if i[0] in sorghum_set:
            rebuilt_fasta_list.append(i[1])
    head(rebuilt_fasta_list)

    output_for_fasta_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\pep_seq_final_1111.fasta", rebuilt_fasta_list)
    output_for_fasta_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\coelorachis_2418_pair_11_11.fasta", new_coelorachis_list)
    # for i in xrange(sorghum_set):





    # for i in
    # seq_return_table = sorted(seq_return_table, key=itemgetter(0))
    # coelorachis_list =
# def set_name():


def seperate_name(astring):
    # astring = astring[1]
    # print astring
    splitted = astring.split(".")
    tail = splitted[-1].split(":")
    splitted = splitted[:-1]
    splitted.extend(tail)
    # exit()
    return splitted


# def load_sorghum():
#     sorghum_file = read_fasta_as_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\Sorghum_3063pair.fasta")


def load_pep_file():
    pep_file = r"Zea_mays.AGPv3.27.pep.all.fa"
    fasta_file = r"D:\2015spring\Bioinfo\Coelorachis\Source data\Zea_mays.AGPv3.27.3325genes_seq7.2.fasta"
    fasta_table = read_fasta_as_table(fasta_file)


    pep_table = read_fasta_as_table(path.join(root_address, pep_file))
    index_list = []
    name_list = [i[0] for i in pep_table]
    name_list = [i.split(":")[0] for i in name_list]
    name_list = [i.split(" ")[0] for i in name_list]
    # name_set = set()
    name_dic = {}
    head(name_list)
    for orderi, i in enumerate(name_list):

        name = i.split("_")[0]
        if name not in name_dic.keys():
            if i[-2:] != "01":continue
            name_dic[name] = orderi
        # else:
        #     name_dic[name].append(orderi)

    name_set = name_dic.keys()
    for orderi, i in enumerate(fasta_table):
        name = i[0]
        name = name.split(".")[2]

        fasta_table[orderi][1] = pep_table[name_dic[name]][1]
        # fasta_table[orderi][0] = pep_table[name_dic[name]][1]
    head(fasta_table)
    print len(fasta_table)
    output_for_fasta_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\Zea_mays.AGPv3.27.3325genes_seq.fasta", fasta_table)







if __name__ == "__main__":
    global root_address
    root_address = r"D:\2015spring\Bioinfo\Coelorachis\Source data"

    # root_address = "D:\\Downloads\\New folder\\"
    # root_address = r"/workdir/bb576/"

    sorghum_coelorachis()
    # load_pep_file()
    # blast_as_ana_pipe()