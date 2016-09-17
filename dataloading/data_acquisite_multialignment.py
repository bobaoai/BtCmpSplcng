from dataloading.data_loading_util import head, sort_list, sort_string, make_fasta_entry_by_list
from dataloading.output_util import *
from dataloading.read_from_gff3 import *
from dataloading.read_from_table import read_table
from dataloading.read_from_fasta import read_fasta_genome_as_dic, read_fasta_as_table
from os import path, walk, rename, makedirs
from Bio import SeqIO
from shutil import copy
__author__ = 'tx'


def fasta_to_nexus(inaddr, outaddr):
    # handle = open("inaddre", "rU")
    # count = SeqIO.convert(inaddr, "fasta", outaddr, "nexus")
    input_handle = open(inaddr, "rU")
    output_handle = open(outaddr, "w")

    sequences = SeqIO.parse(input_handle, "fasta")
    count = SeqIO.write(sequences, output_handle, "nexus")
    output_handle.close()


def convert_fasta_files_in_one_dir():
    output_dir = "multi.nexus"
    output_root = path.join(root_address, output_dir)
    for root, dirs, files in walk("D:\\2015spring\\Bioinfo\\Coelorachis\\multi.afa"):
        for filename in files:
            inaddr = path.join(root, filename)
            print inaddr
            out = filename.split(".afa")
            out = out[0]
            outaddr = path.join(output_root, out + ".nex")
            print outaddr
            fasta_to_nexus(inaddr, outaddr)


def get_seq_from_sorghum(multiseq_dic):
    # input = r"Sorghum_single_scaffold_table_for_seq.p2.7.1.txt"
    # sorghum_table = read_table(path.join(root_address, input))
    # ############
    # load genome

    genome_addr = r"Sorghum_bicolor.Sorbi1.27.dna.genome.fa"
    a_list = read_fasta_as_table(path.join(root_address, genome_addr))
    get_seq_from_sorghum_helper(a_list)
    genome_dic = {}
    for i in a_list:
        genome_dic[i[0]] = i[1]
    # head(a_list)
    ############
    # load scaffold info and gene_length info
    single_scaffold_file = r"Sorghum_single_scaffold_table_for_seq.p3.7.6.txt"
    single_scaffold_table = read_table(path.join(root_address, single_scaffold_file))
    gene_table_with_length = get_gene_from_3425_with_length()
    #############
    # gene_name, scaffold_name, gene_start, gene_end, scaffold_start, scaffold_end, length
    output_fasta_table = []
    for i in xrange(0, len(single_scaffold_table), 2):
        # print
        single_scaffold_table[i].append(gene_table_with_length[single_scaffold_table[i][0]][0])
        k = single_scaffold_table[i]
        first_gene = k[0]
        first_scaffold = k[1]

        init_1, end_1 = get_the_start_end(int(k[2]), int(k[3]), int(k[4]), int(k[5]), int(k[7]))
        # print init, end
        if init_1 <= 0:
            init_1 = 0
        if end_1 > len(genome_dic[k[1]]):
            end_1 = len(genome_dic[k[1]]) - 1

        ############
        single_scaffold_table[i + 1].append(gene_table_with_length[single_scaffold_table[i + 1][0]][0])
        k = single_scaffold_table[i + 1]
        second_gene = k[0]
        second_scaffold = k[1]

        if first_scaffold != second_scaffold:
            print single_scaffold_table[i]
            continue
        init_2, end_2 = get_the_start_end(int(k[2]), int(k[3]), int(k[4]), int(k[5]), int(k[7]))
        # print init, end
        if init_2 <= 0:
            init_2 = 0
        if end_2 > len(genome_dic[k[1]]):
            end_2 = len(genome_dic[k[1]]) - 1

        init = init_1 if init_1 < init_2 else init_2
        end = end_1 if end_1 > end_2 else end_2
        order = first_gene.split(".")
        order = order[0]


        # fasta_name = ":".join([order, first_scaffold, str(init), str(end), "|" + first_gene + "|" + second_gene])
        string = genome_dic[k[1]][init:end + 1]
        multiseq_dic[order].sorghum = ["Sorghum-" + first_scaffold + ":" + str(init) + ":" + str(end), string]

        # output_fasta_table.append([fasta_name, string])

    # output_addr = r"Sorghum_3063pair.fasta"
    # output_for_fasta(output_addr, output_fasta_table)
    return multiseq_dic

    # head(single_scaffold_sorghum_table)
    # head(output_fasta_table)

    # get the name, and truncate it


class maize1():
    def __init__(self, order, tag):
        self.order = order
        self.maize1_gene_name = tag
        self.maize1_gene_start = 0
        self.maize1_gene_end = 0
        self.maize1_chrom = "0"

    def __str__(self):
        return "\t".join([self.order, self.maize1_gene_name,
                          str(self.maize1_gene_start),
                          str(self.maize1_gene_end),
                          str(self.maize1_chrom)])


    def load_chrom_start_end(self, chrom, start, end):
        self.maize1_gene_start = int(start)
        self.maize1_gene_end = int(end)
        self.maize1_chrom = chrom


class multiASeq():
    def __init__(self, order, tag, maize1, maize2):
        temp = tag.split(".")
        self.maize1_gene_name = ".".join(temp[2:])
        self.maize1_gene_start = 0
        self.maize1_gene_end = 0
        self.maize1_chrom = 0

        self.order = order
        self.maize1 = [maize1]
        self.maize2 = [maize2]
        self.sorghum = []
        self.coelorachis = []

    def load_chrom_start_end(self, chrom, start, end):
        self.maize1_gene_start = int(start)
        self.maize1_gene_end = int(end)
        self.maize1_chrom = chrom


def get_seq_from_coelorachis():
    # input = r"Sorghum_single_scaffold_table_for_seq.p2.7.1.txt"
    # sorghum_table = read_table(path.join(root_address, input))
    # ############
    # load genome

    genome_addr = r"Coelorachis.a.lines.fasta"
    a_list = read_fasta_as_table(path.join(root_address, genome_addr))
    # get_seq_from_sorghum_helper(a_list)
    genome_dic = {}
    multiseq_dic = {}
    for i in a_list:
        genome_dic[i[0]] = i[1]
    # exit()
    # for i in genome_dic.keys():
    # print i
    ############
    # load scaffold info and gene_length info
    single_scaffold_file = r"Coelorachis_single_scaffold_table_for_seq.p3.7.6.txt"
    single_scaffold_table = read_table(path.join(root_address, single_scaffold_file))
    gene_table_with_length = get_gene_from_3425_with_length()
    #############
    # gene_name, scaffold_name, gene_start, gene_end, scaffold_start, scaffold_end, length
    output_fasta_table = []
    # Since everything is paired, so we process each pair at one time
    for i in xrange(0, len(single_scaffold_table), 2):
        # print
        single_scaffold_table[i].append(gene_table_with_length[single_scaffold_table[i][0]][0])
        k = single_scaffold_table[i]
        first_gene = k[0]
        first_scaffold = k[1]

        init_1, end_1 = get_the_start_end(int(k[2]), int(k[3]), int(k[4]), int(k[5]), int(k[7]))
        # print init, end
        if init_1 <= 0:
            init_1 = 0
        if end_1 > len(genome_dic[k[1]]):
            end_1 = len(genome_dic[k[1]]) - 1

        ############
        single_scaffold_table[i + 1].append(gene_table_with_length[single_scaffold_table[i + 1][0]][0])
        k = single_scaffold_table[i + 1]
        second_gene = k[0]
        second_scaffold = k[1]

        if first_scaffold != second_scaffold:
            print single_scaffold_table[i]
            continue
        init_2, end_2 = get_the_start_end(int(k[2]), int(k[3]), int(k[4]), int(k[5]), int(k[7]))
        # print init, end
        if init_2 <= 0:
            init_2 = 0
        if end_2 > len(genome_dic[k[1]]):
            end_2 = len(genome_dic[k[1]]) - 1

        init = init_1 if init_1 < init_2 else init_2
        end = end_1 if end_1 > end_2 else end_2
        order = first_gene.split(".")
        order = order[0]
        multiseq_dic[order] = multiASeq(order, first_gene, "Maize1-" + first_gene, "Maize2-" + second_gene)
        multiseq_dic[order].maize1.append(gene_table_with_length[first_gene][1])
        multiseq_dic[order].maize2.append(gene_table_with_length[second_gene][1])

        # fasta_name = ":".join([order, first_scaffold, str(init), str(end), "|" + first_gene + "|" + second_gene])
        string = genome_dic[k[1]][init:end + 1]
        multiseq_dic[order].coelorachis = ["Coelorachis-" + first_scaffold + ":" + str(init) + ":" + str(end), string]

        # output_fasta_table.append([fasta_name, string])

    # output_addr = r"Coelorachis_3063pair.fasta"
    # output_for_fasta(output_addr, output_fasta_table)
    return multiseq_dic

    # head(single_scaffold_sorghum_table)
    # head(output_fasta_table)

    # get the name, and truncate it


def get_gene_from_3425_with_length():
    # return format dic[gene_name] = [sequence length, sequence]
    gene_seq = r"Zea_mays.AGPv3.27.3325.genes_seq_mod_6.26.fasta"
    a_table = read_fasta_as_table(path.join(root_address, gene_seq))
    gene_len_dic = {}
    for i in xrange(len(a_table)):
        gene_len_dic[a_table[i][0]] = [len(a_table[i][1]), a_table[i][1]]

    # head(gene_len_dic, True)
    return gene_len_dic


def get_the_start_end(sca_start, sca_end, gene_start, gene_end, gene_length):
    gene_left = gene_start
    gene_right = gene_length - gene_end
    scaffold_ratio = float(sca_end - sca_start) / (gene_end - gene_start)
    if scaffold_ratio > 20:
        scaffold_left = round(gene_left * 20)
        scaffold_right = round(gene_right * 20)
        scaffold_end = int(sca_end + scaffold_right)
        scaffold_init = int(sca_start - scaffold_left)
        return scaffold_init, scaffold_end

    scaffold_left = round(gene_left * scaffold_ratio)
    scaffold_right = round(gene_right * scaffold_ratio)
    scaffold_end = int(sca_end + scaffold_right)
    scaffold_init = int(sca_start - scaffold_left)
    return scaffold_init, scaffold_end


def get_seq_from_sorghum_helper(a_list):
    for i in xrange(len(a_list)):

        name_string = a_list[i][0]
        name_splited = name_string.split("dna")
        name_splited = name_splited[0][:-1]
        a_list[i][0] = name_splited
        print name_splited


def get_multiaseq_dic():
    multiaseq_dic = get_seq_from_coelorachis()
    multiaseq_dic = get_seq_from_sorghum(multiaseq_dic)
    return multiaseq_dic


def find_coordinates_pipe():
    # multiaseq_dic = get_multiaseq_dic()
    # maize1_name_table = [[multiaseq_dic[i].order, multiaseq_dic[i].maize1_gene_name] for i in multiaseq_dic]
    # head(maize1_name_table)
    # output_for_list(path.join(root_address, "maize1_table"), maize1_name_table)

    # quit()
    gene_name_table = read_table(path.join(root_address, "maize1_table"))
    maize1_dic = {}
    for i in gene_name_table:
        maize1_dic[i[0]] = maize1(i[0], i[1])
    gene_start_end = get_gene_start_end_position_from_gff3(path.join(root_address, "Zea_mays.AGPv3.27.gff3"))

    head(gene_start_end, title="gene_start_end")

    gene_location_dic = {}
    for i in gene_start_end:
        gene_location_dic[i[0]] = i

    for i in maize1_dic.keys():
        chrom, start, end = gene_location_dic[maize1_dic[i].maize1_gene_name][1:4]
        maize1_dic[i].load_chrom_start_end(chrom, start, end)

    v2tov3_table = read_table(path.join(root_address, "v2_v3_map.txt"), head=True)
    head(v2tov3_table, title="v2tov3_table")

    head([str(maize1_dic[i]) for i in maize1_dic.keys()], title="maize_info")
    capture_list = []
    temp_keys = list(maize1_dic.keys())
    deleted = 0
    head(temp_keys)
    for i in v2tov3_table:
        # v2tov3_list.
        temp_capture = CoordinatesCapture(*(i[0:6]))

        dic_keys = list(temp_keys)
        for keys in dic_keys:

            # print temp_capture.v2_chrom, maize1_dic[keys].maize1_chrom
            if temp_capture.v2_chrom != maize1_dic[keys].maize1_chrom: continue
            # print "???"
            start = int(temp_capture.v2_start)
            end = int(temp_capture.v2_end)
            gene_start = maize1_dic[keys].maize1_gene_start
            # print start, end, gene_start
            if start < gene_start < end:
                # print "yes"
                temp_capture.gene_name_list.append(maize1_dic[keys].maize1_gene_name)
                temp_capture.gene_order_list.append(maize1_dic[keys].order)

                temp_keys.remove(keys)
                deleted += 1
        if temp_capture.is_having_gene():
            capture_list.append(temp_capture)
    print deleted
    g = [i.gene_name_list for i in capture_list]
    head(g)
    print "len  gene_name_list", len(g)

    coordinate_file = read_table(path.join(root_address, "maize.paralogous_segments.kNN_10.dat"), head=True)
    coordinate_list = []
    capture_list_copy = list(capture_list)
    total_file = []
    deleted = 0
    for i in coordinate_file:
        temp_coordinate = CoordinateSegment(*i[:3])
        # temp_capture_list = list(capture_list_copy)
        start = temp_coordinate.g1_start
        end = temp_coordinate.g1_end

        for iorder in capture_list_copy:
            if temp_coordinate.g1_chrom != iorder.v1_chrom: continue

            capture_start = iorder.v1_start
            if start < capture_start < end:
                temp_coordinate.capture_gene_list.extend(iorder.gene_name_list)
                temp_coordinate.capture_order_list.extend(iorder.gene_order_list)
                deleted += 1
            else:
                if deleted == 1092:
                    print start, end, capture_start
            # capture_list_copy.remove(iorder)
        # if temp_coordinate.is_having_gene():
        coordinate_list.append(temp_coordinate)
        total_file.append(temp_coordinate.capture_order_list)
    head(total_file)
    print "len total_file", len(total_file)
    print "deleted", deleted
    #################
    # 1216 genes are from coordinate
    #################
    # if len(total_file) == 1476:

    # g = [i.capture_gene_list for i in coordinate_list]
    # head(g)

    # for i in coordinate_list:
    # print(i.capture_gene_list)
    coordinate_dic = {}



    output_dir = "multi.reasigned.nexus"
    a1 = path.join(root_address,output_dir)
    if not path.exists(a1):
        makedirs(a1)
    output_root = path.join(root_address, output_dir)
    output_for_list(path.join(output_root, "check_coordinate"), total_file)
    # for i in xrange(len())
    for i in xrange(len(total_file)):
        coordinate_dic[str(i)] = total_file[i]
        print i, len(total_file[i])
    #     if path.exists(path.join(output_root,str(i))): continue
    #     makedirs(path.join(output_root,str(i)))
    # if not path.exists(path.join(output_root, "outcast")):
    #     makedirs(path.join(output_root, "outcast"))
    #
    #
    # # input_file =
    # for root, dirs, files in walk(path.join(root_address,"multi.nexus")):
    #     for filename in files:
    #         is_outcast = True
    #         input_name = path.join(root, filename)
    #
    #         out = filename.split(".nex")
    #         out = out[0]
    #
    #         # outaddr = path.join(output_root, filename)
    #
    #         for i in coordinate_dic.keys():
    #
    #             if not out in coordinate_dic[i]:continue
    #             if i == str(22):continue
    #             output_temp = path.join(output_root,i)
    #             output_name = path.join(output_temp,filename)
    #             # print input_name, output_name
    #             rename(input_name, output_name)
    #
    #             is_outcast = False
    #         if is_outcast:
    #             output_temp = path.join(output_root,"outcast")
    #             output_name = path.join(output_temp,filename)
    #             # print input_name, output_name
    #
    #             rename(input_name, output_name)



class CoordinateSegment():
    def __init__(self, g1_chrom, g1_start, g1_end):
        self.g1_chrom = g1_chrom
        self.g1_start = int(g1_start)
        self.g1_end = int(g1_end)
        self.capture_gene_list = []
        self.capture_order_list = []

    def is_having_gene(self):
        return self.capture_gene_list != []


class CoordinatesCapture():
    def __init__(self, v1_chrom, v1_start, v1_end, v2_chrom, v2_start, v2_end):
        self.v1_chrom = v1_chrom
        self.v1_start = int(v1_start)
        self.v1_end = int(v1_end)
        self.v2_chrom = v2_chrom
        self.v2_start = int(v2_start)
        self.v2_end = int(v2_end)
        self.gene_name_list = []
        self.gene_order_list = []

    def is_having_gene(self):
        return self.gene_name_list != []


def multialignment_pipe():
    multiaseq_dic = get_multiaseq_dic()
    order_list = list(multiaseq_dic.keys())

    order_list = sort_string(order_list)

    temp_addr = path.join(root_address, "multi")
    output_file = path.join(temp_addr, "order_set")

    order_table = []
    # for i in order_list:
    # order_table.append(r"/programs/muscle/muscle -in multi/"+i+r".fasta -fastaout multi/"+i+r".afa;")
    # head(order_table)
    # output_for_list(output_file, order_table)
    order_table = []
    output_file = path.join(temp_addr, "build_tree")
    for i in order_list:
        order_table.append(r"/programs/FastTree-2.1.7/FastTree -gtr -gamma -nt multi/" + i + ".afa > " + i + ".tre")
    head(order_table)
    output_for_list(output_file, order_table)
    # for i in order_list:
    # output_file = path.join(temp_addr,i+".fasta")
    # # print output_file
    # output_fasta_table = open(output_file,"w")
    # # print make_fasta_entry_by_list(multiaseq_dic[i].maize1)
    #     output_fasta_table.write(make_fasta_entry_by_list(multiaseq_dic[i].maize1))
    #     output_fasta_table.write(make_fasta_entry_by_list(multiaseq_dic[i].maize2))
    #     output_fasta_table.write(make_fasta_entry_by_list(multiaseq_dic[i].coelorachis))
    #     output_fasta_table.write(make_fasta_entry_by_list(multiaseq_dic[i].sorghum))
    #     output_fasta_table.close()


if __name__ == "__main__":
    # read_genome
    global root_address
    root_address = "D:\\2015spring\\Bioinfo\\Coelorachis\\"
    # root_address = "D:\\Downloads\\New folder\\"
    # root_address = r"/workdir/bb576/"
    # global species, version
    # species = "Coelorachis.filtered"
    # species = "Sorghum"
    # version = ".6.22"
    # get_seq_from_sorghum()
    # get_gene_start_end_position_from_gff3()
    multialignment_pipe()
    # convert_fasta_files_in_one_dir()
    find_coordinates_pipe()
