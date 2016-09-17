from operator import itemgetter
from dataloading.data_loading_util import head

__author__ = 'tx'
# import numpy as np
# import matplotlib.pyplot as plt
# from numpy.random import normal
# import plotly.plotly as py
# from plotly.graph_objs import *


def read_table(address, sep="\t", head=False):
    ainput = open(address)

    acc = 0
    alist = []
    if head:
        print ainput.readline()
        pass
        # use dic
    for i in ainput:
        acc += 1
        if i[-1] == "\n":
            i = i[:-1]
        k = i.split(sep)
        alist.append(k)
    print("The file is loaded, and there is " + str(acc) + " lines in total")
    # print "seccessfully loaded the file"
    # for i in xrange(7):
    #     print alist[i]
    return alist






# def save_vector():


def remove_redundance_and_sort(a_list):
    for_remove_same_location = ["\t".join(i) for i in a_list]
    same_location_removed_list = list(set(for_remove_same_location))

    # sort for check

    temp_list = [i.split(".") for i in same_location_removed_list]
    for i in temp_list:
        i[0] = int(i[0])
    temp_list1 = sorted(temp_list)
    for i in temp_list1:
        i[0] = str(i[0])
    # head(temp_list1)
    temp_list = [".".join(i) for i in temp_list1]
    head(temp_list)
    return temp_list


def plot_x(x_list):
    plt.hist(np.asarray(x_list))
    plt.title("Sequence Length")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.show()

def sort_fasta_to_table(addrinput):
    inputfile = open(addrinput)

    geneFastaList = {}
    in_Gene = False
    gene_string = ""
    gene_name = ""
    acc = 0
    for s in inputfile:
        if s.find(">") != -1:
            if gene_name != "":
                # print "hah"
                # if gene_length >= geneLengtThreshold:
                # print "lol"
                if gene_name in geneFastaList.keys():
                    if len(gene_string) >= len(geneFastaList[gene_name]):
                        geneFastaList[gene_name] = gene_string
                else:
                    geneFastaList[gene_name] = gene_string
                if acc < 8:
                    print gene_name, geneFastaList[gene_name]
                    # else:
                    # print "Not Sig"

            gene_name = s[s.find(">") + 1:s.find(" ")]
            if gene_name.find("_") != -1:
                if gene_name.find("_FG") != -1: continue
                gene_name = gene_name[:gene_name.find("_")]
            alist = s.split(":")
            chrom = alist[3]
            gene_start = int(alist[4])
            gene_stop = int(alist[5])
            # print gene_start, gene_stop
            gene_length = gene_stop - gene_start
            gene_string = ""
            acc += 1
        else:
            gene_string += s[:-1]
    # if gene_length >= geneLengtThreshold:
    if gene_name in geneFastaList.keys():
        if len(gene_string) >= len(geneFastaList[gene_name]):
            geneFastaList[gene_name] = gene_string
        else:
            geneFastaList[gene_name] = gene_string

    outputAddr = r"D:\2015spring\Bioinfo\VossiaCoelorachis\FastaTabFile.txt"
    outputFile = open(outputAddr, "w")
    for i in geneFastaList.keys():
        outputFile.write("\t".join([i, geneFastaList[i]]) + "\n")
    outputFile.close()
    return geneFastaList




class tag():
    def __init__(self, **keywords):
        # keywords includes:tag_number, matched_gene, start_position, match_length, gene_startpos
        self._tag_number = keywords['tag_number'] if 'tag_number' in keywords else None
        self._matched_gene = keywords['matched_gene'] if 'matched_gene' in keywords else []
        self._tag_startpos = keywords['start_position'] if 'start_position' in keywords else None
        self._match_length = keywords['match_length'] if 'match_length' in keywords else None
        self._gene_startpos = keywords['gene_startpos'] if 'gene_startpos' in keywords else None
        self._tagpos = ""
        self._tagchrom = ""


    def addtagpos(self, tagpos, tagchrom):
        self.tagpos = tagpos
        self.tagchrom = tagchrom

    @property
    def tagpos(self):
        return self._tagpos

    @tagpos.setter
    def tagpos(self, v):
        self._tagpos = v

    @property
    def tagchrom(self):
        return self._tagchrom

    @tagchrom.setter
    def tagchrom(self, v):
        self._tagchrom = v

    @property
    def tag_number(self):
        return self._tag_number

    @property
    def gene_startpos(self):
        return self._gene_startpos

    @property
    def tag_startpos(self):
        return self._tag_startpos

    @property
    def match_length(self):
        return self._match_length

    @property
    def matched_gene(self):
        return self._matched_gene

    @matched_gene.setter
    def matched_gene(self, string):
        self._matched_gene = string

    def position(self):
        return self.tagchrom + ":" + self.tagpos


class tagGene(tag):
    def __init__(self, **keywords):
        tag.__init__(**keywords)
        self.mappedSeq = {}





def draw_seq_length(aList):
    name_list = [aList[0][-2]]
    len_list = [int(aList[0][-1])]
    for i in aList:
        if i[-2] == name_list[-1]: continue
        if i[-2] in name_list: continue
        name_list.append(i[-2])
        len_list.append(int(i[-1]))
    plt.hist(np.asarray(len_list))
    plt.title("Sequence Length")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.show()























def draw_venn_for_three():
    venn = {}
    venn["Vossia"] = []
    venn["Coelorachis"] = []
    venn["Sorghum"] = []
    for species in ["Vossia", "Coelorachis", "Sorghum"]:
        k = open("D:\\2015spring\\Bioinfo\\VossiaCoelorachis\\single_missed_gene." + species + ".txt")

        for name_string in k:
            gene_name = name_string[:-1]
            venn[species].append(gene_name)

    Vossia_list = list(venn["Vossia"])

    Coelo_list = list(venn["Coelorachis"])
    Sorghum_list = list(venn["Sorghum"])
    head(Vossia_list)
    head(Coelo_list)
    head(Sorghum_list)
    assert len(Vossia_list) == len(list(set(Vossia_list)))
    assert len(Coelo_list) == len(list(set(Coelo_list)))
    assert len(Sorghum_list) == len(list(set(Sorghum_list)))
    a = list(Coelo_list)
    a.extend(list(Vossia_list))
    # head(a)
    Vossia_coel = list(set(a))
    a = list(Sorghum_list)
    a.extend(list(Vossia_list))
    Vossia_Sorghum = list(set(a))
    a = list(Coelo_list)
    a.extend(list(Sorghum_list))
    Sorghum_coel = list(set(a))
    a = list(Sorghum_coel)
    a.extend(list(Vossia_list))
    Vossia_coel_sorghum = list(set(a))

    print len(Vossia_coel), len(Sorghum_coel), len(Vossia_Sorghum), len(Vossia_coel_sorghum)
    all_exist = []
    in_vossia_coel = []
    in_vossia_sorghum = []
    a = list(venn["Vossia"])
    print len(a), len(list(set(a)))
    for i in venn["Vossia"]:
        if i in venn["Coelorachis"] and i in venn["Sorghum"]:
            all_exist.append(i)
        if i in venn["Coelorachis"] and (not i in venn["Sorghum"]):
            in_vossia_coel.append(i)
        if (not i in venn["Coelorachis"]) and i in venn["Sorghum"]:
            in_vossia_sorghum.append(i)
    print len(all_exist), len(in_vossia_sorghum), len(in_vossia_coel)

    all_exist = []
    in_coel_vossia = []
    in_coel_sorghum = []
    a = list(venn["Coelorachis"])
    print len(a), len(list(set(a)))
    for i in venn["Coelorachis"]:
        if i in venn["Vossia"] and i in venn["Sorghum"]:
            all_exist.append(i)
        if i in venn["Vossia"] and (not i in venn["Sorghum"]):
            in_coel_vossia.append(i)
        if (not i in venn["Vossia"]) and i in venn["Sorghum"]:
            in_coel_sorghum.append(i)

    print len(all_exist), len(in_coel_vossia), len(in_coel_sorghum)

    a = list(venn["Sorghum"])
    print len(a), len(list(set(a)))





if __name__ == "__main__":
    # count_aligned_location()
    # pass
    global species
    # for species in ["Vossia", "Coelorachis", "Sorghum"]:
    # print r"D:\2015spring\Bioinfo\VossiaCoelorachis\\"+species+".6.8.filteredAlignments.txt"
    species = "Sorghum"


    # get_single_missed_gene_pipeline()
    # draw_venn_for_three()
