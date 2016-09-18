from _curses import echo

from dataloading.read_from_table import *
from dataloading.read_from_fasta import *


def reads_count_pipeline(ref_data_address, reads_data_address, gene_fasta_file):
    # data from sam
    # data from alignment
    # ref_data_table = read_table(ref_data_address)
    reads_data_table = read_table(reads_data_address)
    print "Number of reads inputted", len(reads_data_table)
    fasta_table = read_fasta_seqio(gene_fasta_file)
    print "Number of transcripts inputted", len(fasta_table)
    # head(fasta_table)

    global gene_dic
    gene_dic = build_gene_name_dic(fasta_table)
    print "dic size", len(list(gene_dic.keys()))

    # load trans's length
    fasta_len = [[i[0], len(i[1])] for i in fasta_table]
    # head(fasta_len)
    get_gene_length(fasta_len)
    # check for trans's length
    lngth_lst = [gene_dic[x].length for x in gene_dic.keys()]
    head(lngth_lst)
    print "Finished length loading"

    # [get_exon_info(x) for x in ref_data_table]

    # load reads_lst
    [clstr_reads_by_gene(x) for x in reads_data_table]
    # check for reads_list
    alst = list(gene_dic.keys())[0:10]
    [head(gene_dic[x].reads_lst) for x in alst]
    print "Finished reads loading"

    [gene_dic[x].count_reads() for x in gene_dic.keys()]
    print "Finished count reads"

    for i in gene_dic.keys():
        gene_dic[i].output_overlap_info()
        # then for some output
    with open(r'/Users/apple/Desktop/BootCampsplicing/counts_output.txt', "w") as text_file:
        for i in gene_dic.keys():

            a_gene = gene_dic[i]
            if a_gene.average_count >= 25: continue
            output_box = ",".join([str(i) for i in a_gene.counts_box_lst])
            output_name = a_gene.name
            output_length = str(a_gene.length)
            output_line = "\t".join([output_name, output_length, output_box])
            text_file.write(output_line + "\n")
    text_file.close()


def get_gene_length(len_lst):
    global gene_dic
    print "start load gene length"
    head(len_lst)
    for i in len_lst:
        gene_dic[i[0]].length = i[1]
        gene_dic[i[0]].set_hit_map()


def build_gene_name_dic(ref_data_table, col=0):
    """
    input: ref_data_table
    *chech the column where gene names are
    :return: set of gene name
    """
    name_list = [i[col] for i in ref_data_table]
    name_set = set(name_list)
    gene_dic = {}
    for i in name_set:
        gene_dic[i] = gene(i)
    return gene_dic


def clstr_reads_by_gene(vector):
    """
    add reads information into gene class
    input:    reads_data_table
    :return:
    """
    global gene_dic
    tran_name, tran_start, tran_end, read_name, length, direction = vector
    mid_point = (float(tran_start) + float(tran_end)) / 2
    if gene_dic[tran_name].length < int(tran_end):
        print "error", tran_name, gene_dic[tran_name].length, tran_start, tran_end
    gene_dic[tran_name].reads_lst.append([int(tran_start), int(tran_end), int(mid_point)])


def get_exon_info(vector):
    """
    add exons information into gene class
    :return:
    """
    global gene_dic
    chrom, read_name, start_pos, end_pos, inverse, gene_name = vector

    pass


class exon:
    def __init__(self, ):
        pass


class gene:
    def __init__(self, g_name):
        self.name = g_name
        self.length = -1
        self.hit_map = None
        self.exon_lst = []
        self.reads_lst = []
        self.counts_box_lst = []
        self.average_count = 0

    def set_hit_map(self):
        if self.length != -1:
            self.hit_map = [0] * self.length
        else:
            print self.length, "already set up?"
            print len(self.hit_map)
            print self.name, "cannot set hit_map"
            # exit()

    def add_read(self):
        pass

    def add_exon(self):
        pass

    def map_reads_on_exon(self):
        pass

    def count_reads(self):
        if self.hit_map is None:
            pass
        for temp_read in self.reads_lst:
            for i in range(temp_read[0], temp_read[1]):
                if i >= len(self.hit_map):
                    print temp_read[1], i, len(self.hit_map), self.name, self.length
                self.hit_map[i] += 1
                # head(self.hit_map)

    def output_overlap_info(self):
        box_num = self.length / 10
        if self.length % 10 != 0:
            box_num += 1
        self.counts_box_lst = [0] * box_num
        for i in xrange(box_num):
            self.counts_box_lst[i] = sum(self.hit_map[i:i + 10]) / 10
        # head(self.counts_box_lst)
        self.average_count = sum(self.counts_box_lst) / len(self.counts_box_lst)
        # exit()


if __name__ == '__main__':
    ref_data_address = ''
    reads_data_address = r'/Users/apple/Desktop/BootCampsplicing/Edahi/Coverage.05.01.bed'
    gene_fasta_file = r"/Users/apple/Desktop/BootCampsplicing/S_pombe_refTrans.fasta"
    reads_count_pipeline(ref_data_address, reads_data_address, gene_fasta_file)
    # unittest.main()
