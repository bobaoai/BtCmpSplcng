from numpy import mean
from blast_result_analysis import blast_file_to_table, filter_identity_ratio
from dataloading.count_aligned_location import GeneCollection, Segment
from dataloading.data_loading_util import head, get_column, sort_list
from dataloading.output_util import output_for_list, output_for_fasta_table
from dataloading.read_from_fasta import read_fasta_as_table
from dataloading.read_from_table import read_table

__author__ = 'tx'


def get_scaffold_set(blast_table):
    return set(get_column(blast_table, 2))


def filter_fasta_file_by_name_list(fasta_file_address, name_list):
    head(name_list)

    fasta_table = read_fasta_as_table(fasta_file_address)
    filtered_fasta_table = []
    for i in fasta_table:
        # if i[0] == "flattened_line_104300":
        # print "we got it"
        # if i[0] in name_list:
        # filtered_fasta_table.append(i)
        # print "gogogo!"
        # name_list.remove(i[0])
        #     else:
        #         print "something wrong!"
        if i[0] in name_list:
            filtered_fasta_table.append(i)
            name_list.remove(i[0])
    # head(filtered_fasta_table)a
    return filtered_fasta_table


class SegmentWithScaffoldInformation(Segment):
    def __init__(self, gene, scaffold_name, gene_start, scaffold_start, length):
        Segment.__init__(self, gene_start, length)
        self.gene_name = gene
        self.scaffold_name = scaffold_name
        self.scaffold_start = scaffold_start
        self.scaffold_end = scaffold_start + length

    def output_table_with_gene_name(self, gene_name):
        return [gene_name, self.scaffold_name, str(self.gene_start), str(self.scaffold_start), str(self.length)]

    def output_blast_string_with_gene_name(self, gene_name=""):
        if gene_name == "":
            return "\t".join([self.gene_name, self.scaffold_name, str(self.gene_start), str(self.scaffold_start),
                              str(self.length)])
        return "\t".join(
            [gene_name, self.scaffold_name, str(self.gene_start), str(self.scaffold_start), str(self.length)])

    def __str__(self):
        print self.scaffold_name, self.gene_name, self.gene_start, self.length, self.repeat_time
        return ""

class GeneVote(GeneCollection):
    """ Merge the best aligned segment and remove others
    """

    def __init__(self, name_string):
        GeneCollection.__init__(self, name_string)
        self.segment_stack = []
        self.aligned_scaffold = {}

    def add_segment_to_stack(self, a_segment):
        # print
        assert isinstance(a_segment, SegmentWithScaffoldInformation)
        self.segment_stack.append(a_segment)

    def merge_aligned_scaffold(self):
        # seg __init__
        # self.gene_start = start
        # self.gene_end = start + length - 1
        # self.length = length
        # self.repeat_time = 1
        # if name:
        # self.name = name
        for a_segment in self.segment_stack:
            if a_segment.scaffold_name not in self.aligned_scaffold.keys():
                self.aligned_scaffold[a_segment.scaffold_name] = a_segment.length
            else:
                self.aligned_scaffold[a_segment.scaffold_name] += a_segment.length

    def merge_collection(self):
        for a_segment in self.segment_stack:
            if not self.collection:
                self.collection.append(a_segment)
            else:
                status, returned_seg_id = self.find_repeat(a_segment)
                if status in [self.TRUE, self.REPLACE]:
                    # print a_segment.scaffold_name
                    # print self.aligned_scaffold[a_segment.scaffold_name]
                    # print self.aligned_scaffold[self.collection[returned_seg_id].scaffold_name]
                    if self.aligned_scaffold[a_segment.scaffold_name] > self.aligned_scaffold[
                        self.collection[returned_seg_id].scaffold_name]:
                        a_segment.repeat_time = self.collection[returned_seg_id].repeat_time
                        self.collection[returned_seg_id] = a_segment
                        # print "lols"
                        # print a_segment.gene_start, a_segment.length, a_segment.repeat_time, a_segment.scaffold_name
                else:
                    self.collection.append(a_segment)

        if self.name =="154.2.GRMZM2G088797":
            for i in self.collection:
                str(i)






    def find_repeat(self, a_segment):
        # modified from the gene_collection
        for stored_segment_id in xrange(len(self.collection)):
            i_segment_gene_start = self.collection[stored_segment_id].gene_start
            i_segment_gene_end = self.collection[stored_segment_id].gene_end
            # aaaaaa
            # bbbbbbb
            if a_segment.gene_start >= i_segment_gene_end or a_segment.gene_end <= i_segment_gene_start:
                continue
                # print "tunnel 1"
                # return self.FALSE

            # aaaaaa
            # bbbbbbbbb
            elif a_segment.gene_start == i_segment_gene_start and a_segment.gene_end == i_segment_gene_end:
                self.collection[stored_segment_id].repeat_time += 1
                # print "tunnel 2"
                return self.TRUE, stored_segment_id
            elif a_segment.gene_start >= i_segment_gene_start and a_segment.gene_end <= i_segment_gene_end:
                self.collection[stored_segment_id].repeat_time += 1
                # print "tunnel 6"
                return self.TRUE, stored_segment_id

            # aaaaaaaaaaaaa
            # bbbbbbbbb
            elif a_segment.gene_start <= i_segment_gene_start and a_segment.gene_end >= i_segment_gene_end:
                self.collection[stored_segment_id].repeat_time += 1
                # print "tunnel 3"
                return self.REPLACE, stored_segment_id

            # aaaaaaaa
            # bbbbbbbbb
            elif a_segment.gene_start <= i_segment_gene_start and a_segment.gene_end <= i_segment_gene_end:
                if self.checkdif(a_segment, self.collection[stored_segment_id]):
                    self.collection[stored_segment_id].repeat_time += 1
                    # print "tunnel 4"
                    return self.TRUE, stored_segment_id
                    # print "tunnel 4"
                    # return self.FALSE

            # aaaaaa
            # bbbbbb
            elif a_segment.gene_start >= i_segment_gene_start and a_segment.gene_end >= i_segment_gene_end:
                if self.checkdif(a_segment, self.collection[stored_segment_id]):
                    self.collection[stored_segment_id].repeat_time += 1
                    # print "tunnel 5"
                    return self.TRUE, stored_segment_id
                    # print "tunnel 5"
                    # return self.FALSE
        return self.FALSE, None

    def output_as_blast_table(self):
        """return table format for blast"""
        output_table = []
        assert isinstance(self.collection[0], SegmentWithScaffoldInformation)
        for a_seg in self.collection:
            output_table.append(a_seg.output_table_with_gene_name(self.name))
        return output_table

    def output_as_blast_string(self):
        output_string_list = []
        assert isinstance(self.collection[0], SegmentWithScaffoldInformation)
        for a_seg in self.collection:
            output_string_list.append(a_seg.output_blast_string_with_gene_name(self.name))
        return output_string_list


def count_aligned_location_by_location_collection(address, output_address, cutoff=00, cutoff_above=1000):
    loaded_list = read_table(address, sep="\t",
                             head=False)
    loaded_list = filter_identity_ratio(loaded_list, 80.0)
    # collect name location_name and length
    a_list = a_list = [[i[1], i[2], int(i[3]), int(i[4]), int(i[5])] for i in loaded_list]
    head(a_list)
    location_collection_dic = {}

    for i in a_list:
        name = i[0]
        # print name
        # start_pos = int(i[1])
        length = int(i[-1])
        # scaffold_name = i[3]
        # scaffold_start = i[4]
        if length < cutoff: continue

        # if length > cutoff_above: continue
        if name not in location_collection_dic.keys():
            location_collection_dic[name] = GeneVote(name)
            location_collection_dic[name].add_segment_to_stack(SegmentWithScaffoldInformation(*i))

        else:
            # print i[1:]
            location_collection_dic[name].add_segment_to_stack(SegmentWithScaffoldInformation(*i))

            # print "add"

    h = 0
    # assert "2.1.GRMZM2G167031" in location_collection_dic.keys(), "noe key found"
    # for i in location_collection_dic["2.1.GRMZM2G167031"].segment_stack:
    #     print str(i)
    # for i in location_collection_dic.keys():
    #     if h < 5:
    #         h += 1
    #         print location_collection_dic[i].segment_stack
    #     else:
    #         break
    # for j in location_collection_dic[i].collection:
    # str(j)
    # h+=1
    repeat_list = []
    for i in location_collection_dic.keys():
        h += 1
        # else: break
        location_collection_dic[i].merge_aligned_scaffold()
        location_collection_dic[i].merge_collection()

        repeat_list.extend(location_collection_dic[i].output_as_blast_string())

    # repeat_list
    # sort the alignment list and output
    a_list = [i.split(".") for i in repeat_list]
    # a_list = [i.split(".") for i in repeat_list]
    a_list = sort_list(a_list)
    # a_list = sorted(repeat_list, key=lambda x: x[0])
    a_list = [".".join(i) for i in a_list]


    # print location_collection_dic["flattened_line_41753"].return_repeat_for_segment()
    output_for_list(output_address,
                    a_list)
    # print mean(repeat_list)


def build_fasta_file_to_filter_redundant_scaffold_pipeline():
    rawBlastTable = blast_file_to_table("/workdir/bb576/Redundacy_removed_Coelorachis.filtered.6.22.txt")
    head(rawBlastTable)
    scaffold_list = list(get_scaffold_set(rawBlastTable))
    inputAddress = "/workdir/bb576/Coelorachis.a.lines.fasta"
    outputAddress = "/workdir/bb576/Coelorachis_filtered.fasta"
    filtered_fasta_table = filter_fasta_file_by_name_list(inputAddress, scaffold_list)
    output_for_fasta_table(outputAddress, filtered_fasta_table)


def check_90identity_ratio(blast_list):
    num_list = [float(i[0]) for i in blast_list]
    posterior = len(num_list)
    prior = sum([i > 88.0 for i in num_list])
    print float(prior) / posterior


def build_voted_fasta_for_multialignment():
    rawBlastTable = blast_file_to_table("D:\\2015spring\\Bioinfo\\Coelorachis\\result_" + species + version + ".txt")
    # filter_identity_ratio(rawBlastTable)
    head(rawBlastTable)
    check_90identity_ratio(rawBlastTable)


if __name__ == "__main__":
    global species, version
    species = "Sorghum.filtered"
    # species = "Sorghum"
    version = ".6.26"
    # species = 'Coelorachis'

    # build_fasta_file_to_filter_redundant_scaffold_pipeline()
    #
    #
    # count_aligned_location_by_location_collection("D:\\2015spring\\Bioinfo\\Coelorachis\\result_"+species+version+".txt", "D:\\2015spring\\Bioinfo\\Coelorachis\\Redundacy_removed_" + species+version + ".txt")
    # # pass

    # build_voted_fasta_for_multialignment()

    # a = [[6438, 107, "a", 2],
    # [6439, 108, "b",3],
    # [6809, 74, "a",4],
    # [6809, 79, "b",5],
    # [6651, 113, "a",6],
    # [6651, 111, "b",7],
    # [6258, 140, "a",8],
    # [6254, 143, "b",9]]
    #
    # a_collection = GeneVote("ab")x
    # for i in a:
    #     a_collection.add_segment_to_stack(SegmentWithScaffoldInformation(i[0], i[1], i[2], i[3]))
    #
    # for i in a_collection.segment_stack:
    #     print i.gene_start, i.length, i.repeat_time, i.scaffold_name
    #
    # a_collection.merge_aligned_scaffold()
    # a_collection.merge_collection()
    # print "break"
    # for i in a_collection.collection:
    #     print i.gene_start, i.length, i.repeat_time, i.scaffold_name
    # print a_collection.output_as_blast_table()
    # aList = blast_file_to_table(rawBlastnFile=r"D:\2015spring\Bioinfo\VossiaCoelorachis\result_" + species + version+".txt")

    # 6.26
    species = "Sorghum.filtered"

    version = ".7.1"
    count_aligned_location_by_location_collection(
        "D:\\2015spring\\Bioinfo\\Coelorachis\\Source data\\Sorghum.6.22.filteredAlignments.txt",
        "D:\\2015spring\\Bioinfo\\Coelorachis\\Redundacy_removed_" + species + version + ".txt")

    species = "Coelorachis.filtered"

    count_aligned_location_by_location_collection(
        "D:\\2015spring\\Bioinfo\\Coelorachis\\Source data\\Coelorachis.filtered.6.22.filteredAlignments.txt",
        "D:\\2015spring\\Bioinfo\\Coelorachis\\Redundacy_removed_" + species + version + ".txt")


