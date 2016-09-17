from numpy import mean
from dataloading.data_loading_util import head
from dataloading.output_util import output_for_list, make_fasta_entry
from dataloading.read_from_table import read_table, remove_redundance_and_sort, plot_x
from read_from_fasta import read_fasta_as_table, read_fasta_as_table1
__author__ = 'tx'


# def count_aligned_location():
# loaded_list = read_table(r"D:\2015spring\Bioinfo\VossiaCoelorachis\Vossia.5.21.filteredAlignments.txt", sep="\t",
# head=False)
#
# a_list = [i[1:3] for i in loaded_list]
#     # head(a_list)
#     list_red_removed_and_sorted = remove_redundance_and_sort(a_list)
#
#     list_red_removed_and_sorted = [i.split("\t") for i in list_red_removed_and_sorted]
#
#     head(list_red_removed_and_sorted)
#     gene_location_count_dic = {}
#     for i in list_red_removed_and_sorted:
#         if i[0] not in gene_location_count_dic.keys():
#             # print i[0]
#             gene_location_count_dic[i[0]] = 1
#         else:
#             gene_location_count_dic[i[0]] += 1
#     get_location_count = []
#     for i in gene_location_count_dic.keys():
#         get_location_count.append(gene_location_count_dic[i])
#     plot_x(get_location_count)
#
#     vec_x = {}
#     head(get_location_count)
#     for i in get_location_count:
#         if str(i) not in vec_x.keys():
#             vec_x[str(i)] = 1
#         else:
#             vec_x[str(i)] += 1
#     veclist = []
#     for key, value in vec_x.iteritems():
#         temp = [int(key), value]
#         veclist.append(temp)
#
#     veclist1 = sorted(veclist)
#     # for i in veclist1:
#     # print i
#     average = weight_average(veclist)
#     out_address = (r"D:\2015spring\Bioinfo\VossiaCoelorachis\Vossia.stat.5.21.txt")
#     standard_output(out_address, get_location_count)
#     print average





# def count_aligned_location_and_add_length():
#     loaded_list = read_table(r"D:\2015spring\Bioinfo\VossiaCoelorachis\Vossia.5.21.filteredAlignments.txt", sep="\t",
#                              head=False)
#
#     # collect name location_name and length
#     a_list = [[i[1], i[2], i[5], i[4]] for i in loaded_list]
#     head(a_list)
#
#     # list_red_removed_and_sorted = remove_redundance_and_sort(a_list)
#     location_collection_dic = {}
#     for i in a_list:
#         name = i[0]
#         location = i[1]
#         length = i[2]
#         position = i[3]
#         if name not in location_collection_dic.keys():
#             location_collection_dic[name] = gene_collection(name)
#         else:
#             if location not in location_collection_dic[name].location_collection.keys():
#                 location_collection_dic[name].location_collection[location] = int(length)
#             else:
#                 location_collection_dic[name].location_collection[location] += int(length)
#
#     location_stat = []
#
#     for i in location_collection_dic.keys():
#         for j in location_collection_dic[i].location_collection.keys():
#             location_stat.append(location_collection_dic[i].location_collection[j])
#
#     # list_red_removed_and_sorted = [i.split("\t") for i in list_red_removed_and_sorted]
#
#     # # head(list_red_removed_and_sorted)
#     # gene_location_count_dic = {}
#     # for i in list_red_removed_and_sorted:
#     # if i[0] not in gene_location_count_dic.keys():
#     #         # print i[0]
#     #         gene_location_count_dic[i[0]] = 1
#     #     else:
#     #         gene_location_count_dic[i[0]] +=1
#     # get_location_count = []
#     # for i in gene_location_count_dic.keys():
#     #     get_location_count.append(gene_location_count_dic[i])
#     # # plot_x(get_location_count)
#     # vec_x = {}
#     # head(get_location_count)
#     # for i in get_location_count:
#     #     if str(i) not in vec_x.keys():
#     #         vec_x[str(i)] = 1
#     #     else:
#     #         vec_x[str(i)] += 1
#     # veclist = []
#     # for key, value in vec_x.iteritems():
#     #     temp = [int(key), value]
#     #     veclist.append(temp)
#     # veclist1 = sorted(veclist)
#     # # for i in veclist1:
#     # #     print i
#     # average = weight_average(veclist)
#     # print average
#     # plot_x(location_stat)
#     # for i in location_stat:
#     #     print i
#     keep_working = True
#     threshold = 50
#     # while(keep_working):
#     #     stat = count_location_per_gene_by_length_cutoff(location_collection_dic,threshold)
#     #     keep_working = True
#     #     temp_mean = mean(stat)
#     #     if temp_mean >3:
#     #         threshold +=5
#     #         continue
#     #     keep_working = False
#     # print threshold, temp_mean
#     stat = count_location_per_gene_by_length_cutoff(location_collection_dic, 150)
#     standard_output(r"D:\2015spring\Bioinfo\VossiaCoelorachis\5.26.txt", stat)
#     print len(stat), mean(stat)


def count_aligned_location_by_location_collection(address, output_address, cutoff=00, cutoff_above=1000):
    loaded_list = read_table(address, sep="\t",
                             head=False)

    # collect name location_name and length
    a_list = [[i[1], i[3], i[5]] for i in loaded_list]
    head(a_list)
    location_collection_dic = {}

    for i in a_list:
        name = i[0]
        start_pos = int(i[1])
        length = int(i[2])
        if length < cutoff: continue

        # if length > cutoff_above: continue
        if name not in location_collection_dic.keys():
            location_collection_dic[name] = GeneCollection(name)
        else:
            location_collection_dic[name].add_segment(Segment(start_pos, length))
            # print "add"
    h = 0
    # for i in location_collection_dic.keys():
    # if h <5: pass
    #         else: break
    #         for j in location_collection_dic[i].collection:
    #             str(j)
    #             h+=1
    repeat_list = []
    for i in location_collection_dic.keys():
        if h < 5: print i
        h += 1
        # else: break

        repeat_list.extend(location_collection_dic[i].return_repeat_for_segment())

    # print location_collection_dic["flattened_line_41753"].return_repeat_for_segment()
    # output_for_list(output_address, repeat_list)

    # used for
    output_for_list(output_address, repeat_list)
    print mean(repeat_list)


def count_location_per_gene_by_length_cutoff(location_collection_dic, cut_length):
    location_stat = []
    for i in location_collection_dic.keys():
        temp_lo = 0
        for j in location_collection_dic[i].location_collection.keys():
            if location_collection_dic[i].location_collection[j] < cut_length: continue
            temp_lo += 1
        location_stat.append(temp_lo)
    return location_stat


def weight_average(a_list):
    total_weight = 0
    total_sum = 0
    for x, y in a_list:
        total_sum += x * y
        total_weight += y
    return float(total_sum) / float(total_weight)


#
# class gene_collection():
# def __init__(self, name, location, length):
#         self.name = name
#         self.location_collection = {}
#         self.location_collection[location]=int(length)


class GeneCollection():
    FALSE = 1
    TRUE = 2
    REPLACE = 2

    def __init__(self, name_string):
        self.collection = []
        self.name = name_string

    def find_repeat(self, a_segment):
        # # print str(a_segment)
        # if self.name == "flattened_line_41753":
        #
        #     print "\n"
        #     print str(a_segment)
        #     for i in self.collection:
        #         print str(i)
        # for i in self.collection:
        #     if i.repeat_time >10:
        #         print self.name
        for stored_segment_id in xrange(len(self.collection)):

            #           aaaaaa
            # bbbbbbb
            if a_segment.gene_start >= self.collection[stored_segment_id].gene_end or a_segment.gene_end <= \
                    self.collection[stored_segment_id].gene_start:
                continue
                # print "tunnel 1"
                # return self.FALSE

            #  aaaaaa
            # bbbbbbbbb
            elif a_segment.gene_start == self.collection[stored_segment_id].gene_start and a_segment.gene_end == \
                    self.collection[stored_segment_id].gene_end:
                self.collection[stored_segment_id].repeat_time += 1
                # print "tunnel 2"
                return self.TRUE
            elif a_segment.gene_start >= self.collection[stored_segment_id].gene_start and a_segment.gene_end <= \
                    self.collection[stored_segment_id].gene_end:
                self.collection[stored_segment_id].repeat_time += 1
                # print "tunnel 6"
                return self.TRUE

            # aaaaaaaaaaaaa
            #   bbbbbbbbb
            elif a_segment.gene_start <= self.collection[stored_segment_id].gene_start and a_segment.gene_end >= \
                    self.collection[stored_segment_id].gene_end:
                a_segment.repeat_time = self.collection[stored_segment_id].repeat_time + 1
                self.collection[stored_segment_id] = a_segment
                # print "tunnel 3"
                return self.REPLACE

            # aaaaaaaa
            #  bbbbbbbbb
            elif a_segment.gene_start <= self.collection[stored_segment_id].gene_start and a_segment.gene_end <= \
                    self.collection[
                        stored_segment_id].gene_end:
                if self.checkdif(a_segment, self.collection[stored_segment_id]):
                    self.collection[stored_segment_id].repeat_time += 1
                    # print "tunnel 4"
                    return self.TRUE
                    # print "tunnel 4"
                    # return self.FALSE

            #    aaaaaa
            # bbbbbb
            elif a_segment.gene_start >= self.collection[stored_segment_id].gene_start and a_segment.gene_end >= \
                    self.collection[stored_segment_id].gene_end:
                if self.checkdif(a_segment, self.collection[stored_segment_id]):
                    self.collection[stored_segment_id].repeat_time += 1
                    # print "tunnel 5"
                    return self.TRUE
                    # print "tunnel 5"
                    # return self.FALSE
        return self.FALSE

    def add_segment(self, a_segment):
        # print
        if not self.collection:
            self.collection.append(a_segment)
        elif self.find_repeat(a_segment) in [self.TRUE, self.REPLACE]:
            pass
        else:
            self.collection.append(a_segment)

    def checkdif(self, a_segment, stored_segment):
        if a_segment.gene_end < stored_segment.gene_end:
            overlap_length = a_segment.gene_end - stored_segment.gene_start
        else:
            overlap_length = stored_segment.gene_end - a_segment.gene_start
        overlap_ratio_for = float(overlap_length)/a_segment.length
        overlap_ratio_for = float(overlap_length)/a_segment.length
        # print a

        # print a_segment.gene_end - stored_segment.gene_end
        # print a
        # ratio = float(a_segment.gene_start - stored_segment.gene_start+a_segment.gene_end - stored_segment.gene_end) / stored_segment.length

        # print a
        # print stored_segment.gene_start
        if self.name =="154.2.GRMZM2G088797" and a_segment.scaffold_name == "3":
            str(stored_segment)
            str(a_segment)
            print overlap_ratio_for
        # print str(a_segment)
        # print str(stored_segment)
        # print abs(ratio)

        if abs(overlap_ratio_for) > 0.6:
            return True
        else:
            return False

    def return_repeat_for_segment(self):
        repeat_list = []
        for i in self.collection:
            repeat_list.append(i.repeat_time)
        return repeat_list


class Segment:
    def __init__(self, gene_start, length):
        self.gene_start = gene_start
        self.gene_end = gene_start + length - 1
        self.length = length
        self.repeat_time = 1


    def __str__(self):
        print self.gene_start, self.length, self.repeat_time

def test_for_6_22_data_source():
    a = read_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\Coelorachis.filtered.6.22.filteredAlignments.txt")
    b = read_fasta_as_table(r"D:\2015spring\Bioinfo\Coelorachis\Source data\coelorachis.lines.fasta")

    scaffold_name_list = [i[2] for i in a]
    scaffold_name_set = set(scaffold_name_list)
    output_file = open(r"D:\2015spring\Bioinfo\Coelorachis\Source data\Coelorachis.filtered.6.25.filteredAlignments.txt", "w")
    for i in b:
        if i[0] in scaffold_name_set:
            output_file.write(make_fasta_entry(i[0], i[1]))
    b = read_fasta_as_table1(r"D:\2015spring\Bioinfo\Coelorachis\Source data\coelorachis.lines.fasta")
    for i in b:
        if i[0] in scaffold_name_set:
            output_file.write(make_fasta_entry(i[0], i[1]))



    # print len(scaffold_name_set)

if __name__ == "__main__":
    # count_aligned_location_and_add_length()
    # count_aligned_location_and_add_length()

    # a = [[6438,107],
    #     [6439, 96],
    #     [6809, 74],
    #     [6809, 75],
    #     [6651,113],
    #     [6651,111],
    #     [6258,140],
    #     [6254,133]]
    #
    # a_collection = location("a")
    # for i in a:
    #
    #     a_collection.add_segment(segment(i[0], i[1]))
    #
    # for i in a_collection.collection:
    #     print i.gene_start,i.length, i.repeat_time
    # global species
    # species = "Sorghum"
    # species = "Vossia"
    # version = ".6.22"
    # # for species in ["Coelorachis", "Sorghum", "CML247"]:
    # output_address = "D:\\2015spring\\Bioinfo\\VossiaCoelorachis\\dis." + species + version + ".txt"
    # count_aligned_location_by_location_collection(
    #     "D:\\2015spring\\Bioinfo\\Coelorachis\\" + species + version + ".filteredAlignments.txt", output_address)

    #11/08
    test_for_6_22_data_source()
