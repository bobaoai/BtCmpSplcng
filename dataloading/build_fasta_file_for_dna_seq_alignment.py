from dataloading.blast_result_analysis import *
from dataloading.data_loading_util import get_column
from dataloading.read_from_gff3 import *
from dataloading.read_from_fasta import *
from dataloading.output_util import *
from dataloading.remove_redundant_scaffold_from_incomplete_genome import filter_fasta_file_by_name_list

__author__ = 'tx'



def get_sequence_from_dic(genome_dic, chrom, start, end):
    return genome_dic[chrom][start:end+1]


def build_gene_nucleotied_seq_fasta_for_blastn():
    # gene_dic = read_fasta_genome_as_dic(
    #     r"/workdir/bb576/Zea_mays.AGPv3.27.dna.genome.fa")
    # gff3_gene_start_end_table = get_gene_start_end_position_from_gff3()
    # gene_table = []
    #
    # for i in gff3_gene_start_end_table:
    #     name = i[0]
    #     chrom_str = i[1]
    #     start = int(i[2])
    #     end = int(i[3])
    #
    #     seq = get_sequence_from_dic(gene_dic, chrom_str, start, end)
    #
    #     gene_table.append([name, seq])
    # head(gene_table )
    # output_for_fasta(r"/workdir/bb576/Zea_mays.AGPv3.27.gene_seq.fa",gene_table)
    rawBlastTable = blast_file_to_table(r"/workdir/bb576/Redundacy_removed_Coelorachis.filtered.6.22.txt")

    [i.insert(0,"") for i in rawBlastTable]
    rawBlastTable = count_gene_amount(rawBlastTable)
    gene_list = list(get_gene_set(rawBlastTable))
    print len(gene_list), "lololol"
    # exit()
    inputAddress = r"/workdir/bb576/Zea_mays.AGPv3.27.3325genes_seq.fasta"
    outputAddress = "/workdir/bb576/Zea_mays.AGPv3.27.3325.genes_seq_mod_6.26.fasta"
    filtered_fasta_table = filter_fasta_file_by_name_list_local(inputAddress, gene_list)
    output_for_fasta(outputAddress, filtered_fasta_table)


def filter_fasta_file_by_name_list_local(fasta_file_address, gene_list_with_order):
    head(gene_list_with_order)

    splited_name_list = [i.split(".") for i in gene_list_with_order]
    gene_list_with_order = [[i[0], i[1], ".".join(i[2:])] for i in splited_name_list]
    head(gene_list_with_order)
    name_list = [".".join(i[2:]) for i in splited_name_list]
    # name_list = list(name_list)
    # print len(name_list), "olololo"
    fasta_table = read_fasta_as_table(fasta_file_address)
    filtered_fasta_table = []
    for i in fasta_table:
        # if i[0] == "flattened_line_104300":
        #     print "we got it"
        #     if i[0] in name_list:
        #         filtered_fasta_table.append(i)
        #         print "gogogo!"
        #         name_list.remove(i[0])
        #     else:
        #         print "something wrong!"
        k = i[0].split(".")
        k = ".".join(k[2:])
        if k in name_list:
            filtered_fasta_table.append([string_in_list(k,gene_list_with_order,2), i[1]])
            name_list.remove(k)
            # print i[0], string_in_list(i[0],gene_list_with_order,2)
    # head(filtered_fasta_table)a
    print len(name_list), "Some gene left?"
    for i in name_list:
        print string_in_list(i,gene_list_with_order,2)
    return filtered_fasta_table

# helper for filter_fasta_local
def string_in_list(string, a_list, column):
    for i in a_list:
        if string == i[column]:
            return ".".join(i)
    assert False, "removed"

def get_scaffold_set(blast_table):
    return set(get_column(blast_table, 2))
def get_gene_set(blast_table):
    return set(get_column(blast_table, 1))

def reform_coelorachis_fasta():
    r"/workdir/bb576/Redundacy_removed_Coelorachis.filtered.6.22.txt"
    rawBlastTable = blast_file_to_table(r"/workdir/bb576/Redundacy_removed_Coelorachis.filtered.6.22.txt")

    [i.insert(0,"") for i in rawBlastTable]
    count_gene_amount(rawBlastTable)
    head(rawBlastTable)
    scaffold_list = list(get_scaffold_set(rawBlastTable))
    inputAddress = "/workdir/bb576/Coelorachis.a.lines.fasta"
    outputAddress = "/workdir/bb576/Coeslorachis_filtered.fasta"
    filtered_fasta_table = filter_fasta_file_by_name_list(inputAddress, scaffold_list)
    output_for_fasta(outputAddress, filtered_fasta_table)

def reform_sorghum_fasta():
    r"/workdir/bb576/Redundacy_removed_Coelorachis.filtered.6.22.txt"
    rawBlastTable = blast_file_to_table(r"/workdir/bb576/Redundacy_removed_Coelorachis.filtered.6.22.txt")

    [i.insert(0,"") for i in rawBlastTable]
    count_gene_amount(rawBlastTable)
    head(rawBlastTable)
    scaffold_list = list(get_scaffold_set(rawBlastTable))
    inputAddress = "/workdir/bb576/Coelorachis.a.lines.fasta"
    outputAddress = "/workdir/bb576/Coeslorachis_filtered.fasta"
    filtered_fasta_table = filter_fasta_file_by_name_list(inputAddress, scaffold_list)
    output_for_fasta(outputAddress, filtered_fasta_table)






if __name__ == "__main__":
    global species, version
    species = "Coelorachis.filtered"
    # species = "Sorghum"
    version = ".6.22"
    # get_gene_start_end_position_from_gff3()
    build_gene_nucleotied_seq_fasta_for_blastn()
    # reform_coelorachis_fasta()
    # read_fasta_as_table(r"/workdir/bb576/Zea_mays.AGPv3.27.3325genes_seq.fasta")