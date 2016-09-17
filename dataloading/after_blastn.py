from dataloading.blast_result_analysis import count_gene_amount
from dataloading.count_aligned_location import GeneCollection
from dataloading.remove_redundant_scaffold_from_incomplete_genome import count_aligned_location_by_location_collection

__author__ = 'tx'
from dataloading.read_from_fasta import *
from dataloading.read_from_table import *
from dataloading.remove_redundant_scaffold_from_incomplete_genome import GeneVote, SegmentWithScaffoldInformation
# build dic for 2 column
from dataloading.output_util import *
from os import path


def build_dic(table_all_string, column_1, column_2, is_int=False):
    a_dic = {}

    if not is_int:
        for i in table_all_string:
            mom = i[column_1]
            son = i[column_2]
            if mom not in a_dic.keys():
                a_dic[mom] = [son]
            else:
                if son in a_dic[mom]:
                    continue
                a_dic[mom].append(son)

    else:
        for i in table_all_string:
            mom = i[column_1]
            son = int(i[column_2])
            if mom not in a_dic.keys():
                a_dic[mom] = [son]
            else:
                if son in a_dic[mom]:
                    continue
                a_dic[mom].append(son)
    head(a_dic, True)
    return a_dic


def build_gene_scaffold_dic(raw_blast_table, gene_column_int=1, scaffold_column_int=2):
    return build_dic(raw_blast_table, gene_column_int, scaffold_column_int)


def filter_by_gene_scaffold_dic(blastn_table, gene_scaffold_dic, gene_column_int=1, scaffold_column_int=2):
    return_table = []
    line_ac = 0
    for i in blastn_table:

        if i[gene_column_int] not in gene_scaffold_dic.keys(): continue

        if i[scaffold_column_int] not in gene_scaffold_dic[i[gene_column_int]]: continue

        return_table.append(i)

    return return_table

    # dic is from cds alignment version


def filter_by_distance_in_same_chrom(a_blasn_table, distance_dic, gene_column_int=1, position_column_int=4,
                                     end_column_int=0, has_end=False):
    return_table = []
    idiot = 0
    for i in a_blasn_table:
        # gene_name = i[gene]
        if i[gene_column_int] not in distance_dic.keys(): continue

        if distance_dic[i[gene_column_int]][0] < distance_dic[i[gene_column_int]][1] - 30000:
            idiot += 1
            print i[gene_column_int]
            print distance_dic[i[gene_column_int]]
            continue
        # gene position = i[position_column_int]
        position = int(i[position_column_int])
        print position, distance_dic[i[gene_column_int]]
        if position < distance_dic[i[gene_column_int]][0] - 5000 or position > distance_dic[
            i[gene_column_int]][1] + 5000: continue
        return_table.append(i)
    print idiot
    return return_table

    # use the early cds alignment version


def build_distance_dic_for_sorghum(raw_blast_table, gene_column_int=1, start_column_int=2):
    start_dic = build_dic(raw_blast_table, gene_column_int, start_column_int, is_int=True)
    return_dic = {}
    for i in start_dic.keys():
        return_dic[i] = []
        for j in start_dic[i]:
            return_dic[i] = build_distance_dic_helper(return_dic[i], j)
            print return_dic[i]
    head(return_dic, True)
    return return_dic


# check location
def build_distance_dic_helper(a_list, position):
    if not a_list:
        return [position, position]
    else:
        if position < a_list[0]:
            a_list[0] = position

        elif position > a_list[1]:
            a_list[1] = position
        return a_list


def filter_for_sorghum():
    cdsBlastTable = read_table(root_address + "Redundacy_removed_Sorghum.filtered.7.1.txt")
    # head(cdsBlastTable)
    ########
    # remove scaffolds which didn't exist in eariler alignment
    a_dic = build_gene_scaffold_dic(cdsBlastTable, gene_column_int=0, scaffold_column_int=1)
    rawBlastTable = read_table(root_address + "result_Sorghum.p3.6.26.txt")
    # print len(rawBlastTable)
    blastn_table_dif_chrom_removed = filter_by_gene_scaffold_dic(rawBlastTable, a_dic)
    # head(blastn_table_dif_chrom_removed)
    # print len(blastn_table_dif_chrom_removed)
    ########
    # remove segment which didn't exist in that position in eariler alignment
    distance_dic = build_distance_dic_for_sorghum(cdsBlastTable, gene_column_int=0, start_column_int=3)
    blastn_table_final_filtered = filter_by_distance_in_same_chrom(blastn_table_dif_chrom_removed, distance_dic)

    head(blastn_table_dif_chrom_removed)
    print len(blastn_table_dif_chrom_removed)
    ####
    # need remove repeat
    ####
    output_for_list(path.join(root_address, "SorghumWithGeneSeq_scaffold_distance_removed.p3.7.6.txt"),
                    blastn_table_final_filtered)

    # get the start and end position on the chromosome
    input_addr = path.join(root_address, "SorghumWithGeneSeq_scaffold_distance_removed.p3.7.6.txt")
    output_addr = path.join(root_address, "SorghumWithGeneSeq_single_scaffold.p3.7.6.txt")
    single_scaffold_capture(input_addr, output_addr)


class GeneMerge(GeneCollection):
    """ Merge the best aligned segment and remove others
    """

    def __init__(self, name_string):
        GeneCollection.__init__(self, name_string)
        self.scaffold_dic = {}
        self.segment_stack = []
        self.start_end_gene = []
        self.start_end_scaffold = []

    def add_segment_to_stack(self, a_segment):
        assert isinstance(a_segment, SegmentWithScaffoldInformation)
        scaffold_name = a_segment.scaffold_name
        if a_segment.scaffold_name not in self.scaffold_dic.keys():
            self.scaffold_dic[scaffold_name] = [a_segment]
        else:
            self.scaffold_dic[scaffold_name].append(a_segment)


    def _add_segment_helper(self, a_list, position):
        if not a_list:
            return [position, position]
        else:
            if position < a_list[0]:
                a_list[0] = position

            elif position > a_list[1]:
                a_list[1] = position
            return a_list

    def is_only_one_scaffold(self):
        return len(self.scaffold_dic.keys()) == 1

    def _merge_scaffold_dic(self):
        # return a table which contain gene_nam, scaffold_name, scaffold_start, scaffold_end
        # merged structure is dic[scaffold_name] = [start, end, length]
        return_dic = {}
        for i in self.scaffold_dic.keys():
            return_dic[i] = self._merge_scaffold_helper(self.scaffold_dic[i])

        return return_dic

    def _merge_scaffold_helper(self, a_segment_stack):
        # self.start_end_gene = []
        # self.start_end_scaffold = []
        scaffold_start_list = [j.scaffold_start for j in a_segment_stack]
        scaffold_end_list = [j.scaffold_end for j in a_segment_stack]
        scaffold_start, scaffold_end = [min(scaffold_start_list), max(scaffold_end_list)]

        start_list = [j.gene_start for j in a_segment_stack]
        end_list = [j.gene_end for j in a_segment_stack]
        gene_start, gene_end = [min(start_list), max(end_list)]

        acc_list = [j.length for j in a_segment_stack]
        acc_length = sum(acc_list)
        # print [scaffold_start, scaffold_end, gene_start, gene_end, acc_length]
        return [scaffold_start, scaffold_end, gene_start, gene_end, acc_length]

    # if they have redundancy or many genes are merged together, it will return []
    def start_end_from_longer_scaffold(self):
        the_merged_scaffold_dic = self._merge_scaffold_dic()
        # head(the_merged_scaffold_dic, True)
        the_key = ""
        the_length = 0
        for i in the_merged_scaffold_dic.keys():

            if the_length < the_merged_scaffold_dic[i][4]:
                the_length = the_merged_scaffold_dic[i][4]
                the_key = i
        a_list = [self.name, the_key]
        a_list.extend(the_merged_scaffold_dic[the_key])
        if a_list[-1] > a_list[-2] - a_list[-3] + 1:
            return []
        # print a_list
        return a_list

    def start_end_from_only_one_scaffold(self):
        if self.is_only_one_scaffold():
            return self.start_end_from_longer_scaffold()
        else:
            assert False, "more than one scaffold"


# check location




def single_scaffold_capture(input, output):
    # input = r"CoelorachisWithGeneSeq_scaffold_redundance_removed.p3.6.26.txt"
    # output = r"CoelorachisWithGeneSeq_single_scaffold.p3.6.26.txt"
    a_file = path.join(root_address, input)
    blastn_table = read_table(a_file)
    if len(blastn_table[0]) == 5:
        [i.insert(0, "") for i in blastn_table]
    # collect name location_name and length

    a_list = [[i[1], i[2], int(i[3]), int(i[4]), int(i[5])] for i in blastn_table]
    head(a_list)
    location_collection_dic = {}

    for i in a_list:
        name = i[0]
        # print name
        # if length > cutoff_above: continue
        if name not in location_collection_dic.keys():
            location_collection_dic[name] = GeneMerge(name)
            location_collection_dic[name].add_segment_to_stack(SegmentWithScaffoldInformation(*i))
        else:
            # print i[1:]
            location_collection_dic[name].add_segment_to_stack(SegmentWithScaffoldInformation(*i))
    # for i in location_collection_dic.keys():
    # location_collection_dic[i].add_segment_to_scaffold_dic()
    gene_count = 0
    single_scaffold_gene_list = []
    head(location_collection_dic, True)
    ######
    # a dictionary with dic[key] = [[geneMerge],[geneMerge]] structure
    gene_pair_dic = {}
    for i in xrange(4509):
        gene_pair_dic[str(i)] = [[], []]
    for i in location_collection_dic.keys():
        # wanna choose gene withsingle scaffold
        # if location_collection_dic[i].is_only_one_scaffold(): continue

        # wanna choose gene with longer scaffold
        # if they have redundancy or many genes are merged together, it will return []
        output_entry = location_collection_dic[i].start_end_from_longer_scaffold()
        # print
        # print output_entry
        if not output_entry:
            continue
        name_split = i.split(".")
        num = (name_split[0])
        order = int(name_split[1]) - 1
        # if
        gene_pair_dic[num][order].append(output_entry)
    paired_count = 0
    # exit()
    output_geneMerge_list = []
    for i in xrange(4509):

        i_str = str(i)
        if gene_pair_dic[i_str][0] == [] or gene_pair_dic[i_str][1] == []:
            # paired_count += 1
            continue
        if len(gene_pair_dic[i_str][0]) > 1 and len(gene_pair_dic[i_str][1]) > 1:
            continue
        if i == 4:
            print gene_pair_dic[i_str][0][0][1]
            print gene_pair_dic[i_str][0][0][0]
        gene_pair_dic[i_str][0], gene_pair_dic[i_str][1] = filter_one_gene_out_from_many_genes(
                gene_pair_dic[i_str][0], gene_pair_dic[i_str][1])

        if gene_pair_dic[i_str][0] == [] or gene_pair_dic[i_str][1] == []:
            # paired_count += 1
            continue

        print gene_pair_dic[i_str][0][0][1]
        if gene_pair_dic[i_str][0][0][1] != gene_pair_dic[i_str][1][0][1]:
            continue

        output_geneMerge_list.extend(gene_pair_dic[i_str][0])
        output_geneMerge_list.extend(gene_pair_dic[i_str][1])

        paired_count += 1
    output_for_list(path.join(root_address, output), output_geneMerge_list, False)




    #########################################################
    # here we have method to output these dic[key] = [[geneMerge],[geneMerge]] structure
    # paired_count = 0
    # output_geneMerge_list = []
    # for i in xrange(4509):
    # i = str(i)
    # if gene_pair_dic[i][0] == [] or gene_pair_dic[i][1] == []:
    # # paired_count += 1
    # continue
    # else:
    # paired_count += 1
    #         output_geneMerge_list.extend(gene_pair_dic[i][0])
    #         output_geneMerge_list.extend(gene_pair_dic[i][1])
    # print paired_count, "paired_count"
    # output_table = []
    # for i in output_geneMerge_list:
    #     for j in i.segment_stack:
    #         output_table.append(j.output_blast_string_with_gene_name(i.name))
    # head(output_table)
    # output_for_list(path.join(root_address, output), output_table)
    ########################################################

    ##########################
    # output for blastn table
    # for i in location_collection_dic.keys():
    # if location_collection_dic[i].is_only_one_scaffold():
    #         for j in location_collection_dic[i].segment_stack:
    #             single_scaffold_gene_list.append(j.output_blast_string_with_gene_name(i))
    #         gene_count += 1
    # output_for_list(path.join(root_address, output), single_scaffold_gene_list)
    ##########################



    # test_list = ["\t"+i for i in single_scaffold_gene_list]

    # count_gene_amount(test_list)

    return single_scaffold_gene_list


def filter_one_gene_out_from_many_genes(gene_list_1, gene_list_2):
    if len(gene_list_1) == 1:
        standard_gene = gene_list_1[0]
        compare_gene_group = gene_list_2
    else:
        standard_gene = gene_list_2[0]
        compare_gene_group = gene_list_1

    overlap_ratio = 0
    return_gene_entry_id = -1
    for i in xrange(len(compare_gene_group)):
        temp_ratio = checkdif(standard_gene, compare_gene_group[i])
        if temp_ratio > overlap_ratio:
            return_gene_entry_id = i
            overlap_ratio = temp_ratio
    if return_gene_entry_id == -1:
        return [], []
    else:
        if len(gene_list_1) == 1:
            return [standard_gene], [compare_gene_group[return_gene_entry_id]]
        else:
            return [compare_gene_group[return_gene_entry_id]], [standard_gene]


# self.name, the_key, scaffold_start, scaffold_end, gene_start, gene_end, acc_lengt
def checkdif(a_entry_1, a_entry_2):
    a_segment_gene_end = a_entry_1[3]
    a_segment_gene_start = a_entry_1[2]
    stored_segment_gene_end = a_entry_2[3]
    stored_segment_gene_start = a_entry_2[2]
    a_segment_length = a_segment_gene_end - a_segment_gene_start
    if a_segment_gene_start < stored_segment_gene_end < a_segment_gene_end:
        pass
    elif a_segment_gene_start < stored_segment_gene_start < a_segment_gene_end:
        pass
    else:
        return 0

    if a_segment_gene_end < stored_segment_gene_end:
        overlap_length = a_segment_gene_end - stored_segment_gene_start
    else:
        overlap_length = stored_segment_gene_end - a_segment_gene_start
    overlap_ratio_for = float(overlap_length) / a_segment_length
    assert overlap_ratio_for > 0
    return overlap_ratio_for
    # print a

    # print a_segment.gene_end - stored_segment.gene_end
    # print a
    # ratio = float(a_segment.gene_start - stored_segment.gene_start+a_segment.gene_end - stored_segment.gene_end) / stored_segment.length

    # print a
    # print stored_segment_gene_start
    # if self_name == "154.2.GRMZM2G088797" and a_segment.scaffold_name == "3":
    # str(stored_segment)
    #     str(a_segment)
    #     print overlap_ratio_for
    # print str(a_segment)
    # print str(stored_segment)
    # print abs(ratio)


def filter_for_coelorachis():
    cdsBlastTable = read_table(path.join(root_address, "Redundacy_removed_Coelorachis.filtered.7.1.txt"))
    head(cdsBlastTable)


    ########
    # remove scaffolds which didn't exist in eariler alignment
    a_dic = build_gene_scaffold_dic(cdsBlastTable, gene_column_int=0, scaffold_column_int=1)
    rawBlastTable = read_table(path.join(root_address, "result_Coelorachis.p3.6.26.txt"))
    # head(rawBlastTable)
    blastn_table_scaffold_removed = filter_by_gene_scaffold_dic(rawBlastTable, a_dic)
    # head(blastn_table_scaffold_removed)

    output_for_list(path.join(root_address, "Coelorachis_unrelaventScaffold_removed.p3.7.6.txt"),
                    blastn_table_scaffold_removed)
    ########
    # removed overlapped part
    count_aligned_location_by_location_collection(
        path.join(root_address, "Coelorachis_unrelaventScaffold_removed.p3.7.6.txt"),
        path.join(root_address,
                  "CoelorachisWithGeneSeq_scaffold_redundance_removed.p3.7.6.txt"))
    #########
    # find the longest part of alignment
    input_addr = r"CoelorachisWithGeneSeq_scaffold_redundance_removed.p3.7.6.txt"
    output_addr = r"CoelorachisWithGeneSeq_single_scaffold.p3.7.6.txt"
    single_scaffold_capture(input_addr, output_addr)


def get_common_gene_from_sorghum_and_coelorachis():
    input_file1 = r"CoelorachisWithGeneSeq_single_scaffold.p3.7.6.txt"
    input_file2 = r"SorghumWithGeneSeq_single_scaffold.p3.7.6.txt"

    output_file1 = r"Coelorachis_single_scaffold_table_for_seq.p3.7.6.txt"
    output_file2 = r"Sorghum_single_scaffold_table_for_seq.p3.7.6.txt"

    coelorahics_table = read_table(path.join(root_address, input_file1))
    sorghum_table = read_table(path.join(root_address, input_file2))
    # coelorahics_num_order_name = [i[0] for i in coelorahics_table]
    # sorghum_num_order_name = [i[0] for i in sorghum_table]
    #
    # coelorahics_num = [[i[0], i[1]] for i in coelorahics_num_order_name]
    # sorghum_num = [[i[0], i[1]] for i in sorghum_num_order_name]

    # coelorahics_num = list(set(coelorahics_num))
    # sorghum_num = list(set(sorghum_num))
    coelorachis_list = []
    sorghum_list = []
    for i in coelorahics_table:

        for j in sorghum_table:
            if i[0] == j[0]:
                coelorachis_list.append(i)
                sorghum_list.append(j)

    coelorachis_output = count_gene_amount(coelorachis_list, with_identity_ratio=False)
    # head(coelorachis_output)
    sorghum_output = count_gene_amount(sorghum_list, with_identity_ratio=False)
    # head(sorghum_output)

    output_for_list(path.join(root_address, output_file1), coelorachis_output)
    output_for_list(path.join(root_address, output_file2), sorghum_output)


def after_blastn_pipeline():
    # we need two files as input
    filter_for_coelorachis()
    filter_for_sorghum()
    get_common_gene_from_sorghum_and_coelorachis()


if __name__ == "__main__":
    global root_address
    root_address = "D:\\2015spring\\Bioinfo\\Coelorachis\\"
    # root_address = "D:\\Downloads\\New folder\\"
    # root_address = r"/workdir/bb576/"

    after_blastn_pipeline()



    # a = [[],[]]
    # print a
