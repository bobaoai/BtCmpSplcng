from dataloading.data_loading_util import head

__author__ = 'tx'

from Bio import SeqIO

def read_fasta_title(address):
    file = open(address)
    name_list = []
    for i in file:
        if i.find(">") != -1:
            k = i[:i.find("\n")]
            name_list.append(k)
    return name_list



def read_fasta_seqio(address):
    fasta_sequences = SeqIO.parse(open(address),'fasta')
    return_list = []
    for seq_record in fasta_sequences:
        # print(seq_record.id)
        # print(repr(seq_record.seq))
        # print(len(seq_record))
        return_list.append([seq_record.id, seq_record.seq])
    return return_list


def read_fasta_as_dic(address):
    a_list = read_fasta_as_table(address)
    return_dic = {}
    for i in a_list:
        return_dic[i[0]] = i[1]
    return return_dic


def read_fasta_as_table(address):
    opened_file = open(address)
    # prepare attri start
    in_line = False
    fasta_table = []
    string_acc = 0
    seq_str = ""
    seq_name = ""
    acc = 0
    # prepare attri end

    for order, i in enumerate(opened_file):

        i = i[:-1]
        if order % 1000000==0:
            print order
        #     break
        # if acc > 100:
        #     break
        if i.find(">") != -1:

            acc += 1
            if seq_name != "":

                fasta_table.append([seq_name, seq_str])

            seq_name = i[1:]
            # print "Start fasta: ", seq_name
            seq_str = ""
        else:
            seq_str += i
    fasta_table.append([seq_name, seq_str])

    # head(fasta_table)
    print "Successfully load fasta file, there are " + str(len(fasta_table)) + " lines"
    return fasta_table

def read_fasta_as_table1(address):
    opened_file = open(address)
    # prepare attri start
    in_line = False
    fasta_table = []
    string_acc = 0
    seq_str = ""
    seq_name = ""
    acc = 0
    # prepare attri end

    for order, i in enumerate(opened_file):

        i = i[:-1]
        if order < 500127:
                continue
        # if acc > 100:
        #     break
        if i.find(">") != -1:
            # print order
            # exit()
            acc += 1
            if seq_name != "":

                fasta_table.append([seq_name, seq_str])

            seq_name = i[1:]
            # print "Start fasta: ", seq_name
            seq_str = ""
        else:
            seq_str += i
    fasta_table.append([seq_name, seq_str])

    # head(fasta_table)
    print "Successfully load fasta file, there are " + str(len(fasta_table)) + " lines"
    return fasta_table


def read_fasta_genome_as_dic(genome_file):
    # _read_genome: loading the genome seq from parameter: genome_file
    #   into parameter list: chrom which is a list of class: chrom_cds
    genome_handle = open(genome_file)
    chrom = {}
    LIST = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'Mt', 'UNKNOWN', 'Pt']
    for i in LIST:
        chrom[i] = ""
    seq_name = ""
    chrom_name =""
    seq_str = ""
    end_position_int = 0
    for a in genome_handle:
        a = a[:-1]
        if a[0] == '>':
            # if a[2]!="1":
            #    b=1
            #    break

            if chrom_name != "":
                chrom[chrom_name] = seq_str
                print len(seq_str)
                assert len(seq_str) == end_position_int
            chrom_name =""
            seq_name = a[1:a.find(' ')]
            seq_str = ""
            if seq_name in LIST:
                chrom_name = seq_name
                a_splited_str = a.split("AGPv3:"+seq_name+":1:")
                end_position_int = int(a_splited_str[1].split(":")[0])
                # print a
                print "start seq: " + seq_name
                print "Length is " + str(end_position_int)
        else:
            if chrom_name in LIST:
                seq_str += a
    if seq_name in LIST:
        chrom[seq_name] = seq_str
        assert len(seq_str) == end_position_int
    # for i in chrom:
    #    j=''
    #    i.seq=''.join(i.seq)
    print "Finish_loading_genome"
    return chrom

# if __name__ == "__main__":
# address = r"/workdir/bb576/CML_500bp_S1_L001_R1_001.fastq"
# out = r"/workdir/bb576/CML_500bp_S1_L001_R1_001.parse.fastq"
# get_line()
# address =
