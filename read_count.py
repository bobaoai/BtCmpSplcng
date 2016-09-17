import unittest

from dataloading import read_from_table


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)


def reads_count_pipeline(ref_data_address, reads_data_address):
    # data from sam
    # data from alignment
    ref_data_table = read_from_table(ref_data_address)
    reads_data_table = read_from_table(reads_data_address)
    global gene_dic
    gene_dic = build_gene_name_dic(ref_data_table)
    [get_exon_info(x) for x in ref_data_table]
    [clstr_reads_by_gene(x) for x in reads_data_table]
    [gene_dic[x].count_reads() for x in gene_dic.keys()]

    for i in gene_dic.keys():
        gene_dic[x].output_overlap_info()
        #then for some output



def build_gene_name_dic(ref_data_table, col=6):
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


def clstr_reads_by_gene():
    """
    add reads information into gene class
    input:    reads_data_table
    :return:
    """
    global gene_dic

    pass


def get_exon_info():
    """
    add exons information into gene class
    :return:
    """
    global gene_dic
    pass


class exon:
    def __init__(self, ):
        pass


class gene:
    def __init__(self, g_name):
        self.name = g_name
        self.exon_lst = []
        self.reads_lst = []

    def add_read(self):
        pass

    def add_exon(self):
        pass

    def map_reads_on_exon(self):
        pass

    def count_reads(self):
        pass
    def output_overlap_info(self):
        pass


if __name__ == '__main__':
    ref_data_address = ''
    reads_data_address = ''
    reads_count_pipeline(ref_data_address, reads_data_address)
    # unittest.main()
