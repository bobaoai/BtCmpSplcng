from operator import itemgetter

__author__ = 'tx'


def head(f, is_dic=False, title=""):
    if title != "":
        print title
    print "HEAD "
    k = 0
    if is_dic:
        print "Head for dic"
        line_ac = 0
        for i in f.keys():
            if f[i] is list:
                print i + ":"
                head(f[i])
            else:
                print i
                print f[i]
            line_ac += 1
            if line_ac == 10:
                return
    else:
        for i in f:
            print i
            if k > 10:
                break
            k += 1


def make_fasta_entry(name, seq):
    return ">" + name + "\n" + seq + "\n"
def make_fasta_entry_by_list(a_list):
    name = a_list[0]
    seq = a_list[1]
    return ">" + name + "\n" + seq + "\n"


class list_iter_util():
    def __init__(self):
        self.acc = 0

    def add_one(self):
        self.acc += 1

    def larger_than(self, time):
        if time < self.acc:
            return True
        else:
            return False


def get_column(aList, s):
    if type(s) == int:
        return [i[s] for i in aList]
    if type(s) is list:
        new_list = []
        for i in aList:
            temp = []
            for j in s:
                temp.append(i[j])
            new_list.append(temp)
        return new_list


def sort_list(temp_list):
    for i in temp_list:
        i[0] = int(i[0])
    temp_list1 = sorted(temp_list, key=itemgetter(0))
    for i in temp_list1:
        i[0] = str(i[0])
    return temp_list1

def sort_string(temp_list):
    string_list = []
    for i in temp_list:
        string_list.append(int(i))
    sortted_list = sorted(string_list)
    return_list = []
    for i in sortted_list:
        return_list.append(str(i))
    return return_list