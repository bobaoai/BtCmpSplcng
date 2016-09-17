from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from dataloading.output_util import output_for_list, output_for_fasta_dic

__author__ = 'tx'
from os import path, walk
from Bio import SeqIO
import re
from Bio.Phylo.Consensus import *
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO

root_address = "D:\\2015spring\\Bioinfo\\Coelorachis\\"
ATCG = re.compile("--*")


def fasta_to_format(inaddr="", outaddr="", format="aln"):
    # def fasta_to_format(format="aln"):
    # handle = open("inaddre", "rU")
    # count = SeqIO.convert(inaddr, "fasta", outaddr, "nexus")
    # inaddr = "D:\\2015spring\\Bioinfo\\Coelorachis\\merged_multifasta500.afa"
    # outaddr = "D:\\2015spring\\Bioinfo\\Coelorachis\\merged_multifasta500.nex"

    input_handle = open(inaddr, "rU")
    # output_handle = open(outaddr, "w")
    # output_handle.close()

    sequences = SeqIO.parse(input_handle, "fasta")
    # for i in sequences:
    # print i.id
    # count = 0

    a = list(sequences)

    for j in a:
        # sequence = j.seq

        sequence = str(j.seq)
        # print sequence
        if sequence.find("N") != -1:
            # count += 1
            have_N = True
            return
    for k in a:
        if k.id == "Coelorachis":
            sequence = k.seq
            break
    precount = 0
    while sequence[0] == "-":
        sequence = sequence[1:]
        precount += 1
    postcount = 0
    while sequence[-1] == "-":
        sequence = sequence[0:-1]
        postcount -= 1
    if postcount == 0:

        for i in a:
            i.seq = i.seq[precount:]
    else:
        for i in a:
            i.seq = i.seq[precount: postcount]
    # if have_N:
    # count2 += 1
    #     continue
    output_handle = open(outaddr, "w")

    count = SeqIO.write(a, output_handle, "phylip")
    output_handle.close()
    # try:
    #     multi_tree(outaddr)
    # except ZeroDivisionError:
    #     pass


def lola():
    outaddr = "D:\\2015spring\\Bioinfo\\Coelorachis\\merged_multifasta.phy"
    aln = AlignIO.read(outaddr, 'phylip')
    calculator = DistanceCalculator('blastn')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    tree.root_with_outgroup("Sorghum")
    print str(tree.format("newick"))


def multi_tree(outaddr):
    # outaddr = "D:\\2015spring\\Bioinfo\\Coelorachis\\merged_multifasta.phy"
    aln = AlignIO.read(outaddr, 'phylip')
    calculator = DistanceCalculator('blastn')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    tree.root_with_outgroup("Sorghum")
    print str(tree.format("newick"))


def get_file_list_and_out_list(addr, output_root, inputformat="", outputformat=""):
    return_list = []

    for root, dirs, files in walk(addr):
        for filename in files:
            inaddr = path.join(root, filename)
            # print inaddr
            out = filename
            if inputformat != "":
                out = filename.split(".afa")
                out = out[0]
            if outputformat != "":
                outname = out + outputformat
            else:
                outname = filename
            outaddr = path.join(output_root, outname)
            # print outaddr
            return_list.append([inaddr, outaddr])
    return return_list


def merge_msa_fasta_file():
    addr = "D:\\2015spring\\Bioinfo\\Coelorachis\\multiafa.namechanged"
    out = "D:\\2015spring\\Bioinfo\\Coelorachis\\merged_multifasta500.afa"
    name_file = get_file_list_and_out_list(addr, out)
    name_file = [i[0] for i in name_file]

    seq_dic = {"Maize1": "", "Coelorachis": "", "Maize2": "", "Sorghum": ""}
    # record = SeqIO.parse(name_file[0], "fasta")
    count = 0
    count2 = 0
    # have_N = False
    gene_location_list = [0]
    for order, i in enumerate(name_file):
        if order > 500:
            break
        have_N = False
        # if order ==10:
        # break
        print order
        record = SeqIO.parse(i, "fasta")
        a = list(record)
        for j in a:
            # sequence = j.seq

            sequence = str(j.seq)
            # print sequence
            if sequence.find("N") != -1:
                count += 1
                have_N = True
                break
        if have_N:
            count2 += 1
            continue
        for j in a:
            # print j.id
            # print j.seq
            seq_dic[j.id] += str(j.seq)
        gene_location_list.append(len(seq_dic["Sorghum"]))
        # print len(seq_dic[j.id])
    # for i in seq_dic.keys():
    # print seq_dic[i]
    print count
    print count2
    output_for_fasta_dic(out, seq_dic)
    # print len(seq_dic["Sorghum"])
    # output_nex_info = open("D:\\2015spring\\Bioinfo\\Coelorachis\\nexus_info.txt", "w")
    # output_nex_info.write("begin mrbayes;\nset autoclose=yes nowarn=yes;\noutgroup 4;\n")
    # gene_name_list = []
    # for i in xrange(len(gene_location_list) - 1):
    # output_nex_info.write("CHARSET gene" + str(i) + " = " + str(gene_location_list[i]+1) + "-" + str(
    #         gene_location_list[i+1]) + ";\n")
    #     gene_name_list.append("gene" + str(i))
    #
    # output_nex_info.write("partition Genes = " + str(len(gene_location_list)-1) + ": ")
    # output_nex_info.write(",".join(gene_name_list))
    # output_nex_info.write(";\n")
    # output_nex_info.write("set partition = Genes;\ntaxset Maize2 = Maize2;\ntaxset Coelorachis = Coelorachis;\ntaxset Maize1 = Maize1;\ntaxset Sorghum = Sorghum;")
    # output_nex_info.write("lset nst=6 rates=invgamma;\nprset thetapr=invgamma(3,0.003) GeneMuPr=uniform(0.5,1.5) best=1;\n"
    #                       "unlink topology=(all) brlens=(all)  genemu=(all) shape=(all) Statefreq=(all)tratio=(all) revmat=(all) pinvar=(all);\n"
    #                       "mcmc ngen=20000000 samplefreq=2000 nrun=4 nchain=4;\nquit;\nend;\n")


def convert_fasta_files_in_one_dir():
    output_dir = "multiphylip.namechanged"
    output_root = path.join(root_address, output_dir)
    for root, dirs, files in walk("D:\\2015spring\\Bioinfo\\Coelorachis\\multiafa.namechanged"):
        for filename in files:
            inaddr = path.join(root, filename)
            # print inaddr
            out = filename.split(".afa")
            out = out[0]
            outaddr = path.join(output_root, out + ".phylip")
            # print outaddr
            fasta_to_format(inaddr, outaddr)


maize2 = re.compile("Maize2-[^\n]*")
Sorghum = re.compile("Sorghum-[^\n]*")
maize1 = re.compile("Maize1-[^\n]*")
Coelorachis = re.compile("Coelorachis-[^\n]*")


def _change_leaf_name_helper(inaddr, outaddr):
    file = list(open(inaddr))

    output = open(outaddr, "w")
    for tree in file:
        tree = re.sub(maize1, "Maize1", tree)
        tree = re.sub(Sorghum, "Sorghum", tree)
        tree = re.sub(maize2, "Maize2", tree)
        tree = re.sub(Coelorachis, "Coelorachis", tree)
        output.write(tree)
    output.close()


def change_leaf_name():
    output_root = "D:\\2015spring\\Bioinfo\\Coelorachis\\multiafa.namechanged"
    for root, dir, files in walk("D:\\2015spring\\Bioinfo\\Coelorachis\\multi.afa"):
        for filename in files:
            inaddr = path.join(root, filename)
            print inaddr
            # out = filename.split(".tre")
            # out = out[0]
            outaddr = path.join(output_root, filename)
            print outaddr
            _change_leaf_name_helper(inaddr, outaddr)


def change_root_name():
    output_root = "D:\\2015spring\\Bioinfo\\Coelorachis\\multitre.namechanged.rerooted"
    newick_tree = []
    trees_output_addr = "D:\\2015spring\\Bioinfo\\Coelorachis\\multitre.newicktrees"
    for root, dir, files in walk("D:\\2015spring\\Bioinfo\\Coelorachis\\multitre.namechanged"):
        for filename in files:
            inaddr = path.join(root, filename)
            print inaddr
            # out = filename.split(".tre")
            # out = out[0]
            outaddr = path.join(output_root, filename)
            print outaddr
            newick_tree.extend(_change_root_helper(inaddr, outaddr))
    output_for_list(trees_output_addr, newick_tree)


def _change_root_helper(inaddr, outaddr):
    trees = Phylo.parse(inaddr, "newick")
    output_file = open(outaddr, "w")
    newick_tree = []
    for tree in trees:
        tree.root_with_outgroup("Sorghum")
        newick_tree.append(str(tree.format("newick")))
        output_file.write(str(tree.format("newick")))
    output_file.close()

    return newick_tree
    # quit()


def consensus_tree():
    trees_input_addr = "D:\\2015spring\\Bioinfo\\Coelorachis\\merged_multitree_newick.txt"
    trees = list(Phylo.parse(trees_input_addr, "newick"))
    print len(trees)
    # for tree in trees:
    # print str(tree.format("newick"))
    strict_tree = strict_consensus(trees)
    print str(strict_tree.format("newick"))
    Phylo.draw_ascii(strict_tree)

    majority_tree = majority_consensus(trees, 0.2)
    print str(majority_tree.format("newick"))
    Phylo.draw_ascii(majority_tree)


    adam_tree = adam_consensus(trees)
    print str(adam_tree.format("newick"))
    Phylo.draw_ascii(adam_tree)


if __name__ == "__main__":
    # read_genome
    # pass
    # change_leaf_name()
    # change_root_name()
    # consensus_tree()
    convert_fasta_files_in_one_dir()
    # fasta_to_format()
    # merge_msa_fasta_file()
    # lola()