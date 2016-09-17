from dataloading.data_loading_util import head
from dataloading.read_from_table import read_table

__author__ = 'tx'

def read_gff3_file_to_table(address):
    gff3_table=read_table(address, sep ="\t", head=False)
    remove_preface(gff3_table)
    # head(gff3_table)
    gff3_table = remove_for_only_nucle_genome_info(gff3_table)
    return gff3_table
# def
def remove_preface(a_table):
    is_Preface = True
    line = 0
    while is_Preface:
        if a_table[line][0][0:2] == "##":
            a_table.pop(line)
        else:
            is_Preface = False

def remove_for_only_nucle_genome_info(a_table):
    return_list =[]
    chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Mt', 'UNKNOWN', 'Pt']

    print(chrom)
    for i in a_table:
        if not i[0] in chrom:continue
        return_list.append(i)
    return return_list


def get_gene_entry(a_gff_entry):
    gene_mess = a_gff_entry[8]
    gene_mess_splited=gene_mess.split("ID=gene:")
    gene_mess_splited=gene_mess_splited[1].split(";")
    gene_name = gene_mess_splited[0]
    i = a_gff_entry
    return [gene_name, i[0] ,i[3], i[4], i[6]]


def get_gene_start_end_position_from_gff3(addr=r"/workdir/bb576/Zea_mays.AGPv3.30.gff3"):
    gff3_table = read_gff3_file_to_table(addr)
    gene_gff3_table = gff3_gene_entry_table(gff3_table)
    # head(gene_gff3_table)
    info_table = [get_gene_entry(i) for i in gene_gff3_table]
    head(info_table)
    return info_table


def gff3_gene_entry_table(a_gff3_table):
    return_list = []
    for i in a_gff3_table:
        if i[2] == "gene":
            return_list.append(i)
    return return_list



if __name__ == "__main__":
    global species, version
    species = "Coelorachis.filtered"
    # species = "Sorghum"
    version = ".6.22"
    read_gff3_file_to_table(r"D:\2015spring\Bioinfo\Coelorachis\Zea_mays.AGPv3.30.gff3")
    # species = 'Coelorac
