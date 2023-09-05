import argparse
import time
from collections import deque
from collections import defaultdict
from collections import Counter

import Levenshtein
import numpy as np


# 目的基于minizer识别不同长度和错配允许个数下的dup
# 对于一个整体的fastq文件，应该有一个整体的长度值/均值，同时错配数+1是windows数目，windows长度-错配数是k值。举个例子，150长度允许4处错误就是，30w,26k,3处错误就是37，37，38，38,不同区域分别37-3和38-3，120长度3处错误，30w-27k



def find_dup_base_distance(all_read,dis=4,num_for_diff=8):

    if len(all_read)<10:
        all_sci_dir = defaultdict(list)
        for i in range(len(all_read)):
            for j in range(i+1,len(all_read)):
                if Levenshtein.distance(all_read[i], all_read[j]) < dis:
                    all_sci_dir[str(i)].append(j)
        return all_sci_dir

    all_del_pos_list = np.triu(np.ones((len(all_read), len(all_read)), dtype=bool), k=1)
    all_sci_dir = {}

    for read_num in range(num_for_diff):
        temp_line = np.array([Levenshtein.distance(all_read[read_num], s2) for s2 in all_read], dtype=np.int16)
        temp_condition = np.abs(np.subtract.outer(temp_line, temp_line)) < dis
        all_del_pos_list &= temp_condition
        temp_aa = np.where(temp_line < dis)[0]
        if len(temp_aa) != 0:
            all_sci_dir[str(read_num)] = temp_aa

    for read_num in range(num_for_diff, len(all_read)):
        temp_index = np.where(all_del_pos_list[read_num] == 1)[0].astype(np.uint32)
        temp_line = np.zeros(len(temp_index), dtype=np.int16)
        for i, my_index in enumerate(temp_index):
            temp_line[i] = Levenshtein.distance(all_read[read_num], all_read[my_index])

        temp_aa = np.where(temp_line < dis)[0]
        if len(temp_aa) != 0:
            all_sci_dir[str(read_num)] = temp_index[temp_aa]

    return all_sci_dir

# 对于读取的fastq文件，如何保存，保存name.split()[1],seq,quality,保存到什么数据结构   ，肯定是字典，key是kmer，value是name，其他都是中间内容

# 我的建议是别差这点内存，额外再来两个字典，key是name，value是seq,  不好，还是字符串空格吧，元组后面再试

def is_substring_with_one_error(main_string1,main_string2, substring1="ATGTGTATAAGA",substring2="GTGTATAAGAGA",substring3="ATATATAAAAAA"):
    result=False
    start1 = main_string1.find(substring1)
    start2 = main_string2.find(substring1)
    barcode1=main_string1[start1-10:start1-2]
    barcode2=main_string2[start2-10:start2-2]
    if barcode1.replace("C","T")==barcode2.replace("C","T"):
        result=True

    return result


def get_barcode(main_string, substring):
    start = main_string.find(substring)
    if start != -1:

        return start, start + len(substring)


def judge(str1, str2):
    order = {'G': 0, 'A': 1, 'T': 2, 'N':5}
    for char1, char2 in zip(str1, str2):
        if order[char1] < order[char2]:
            return False
        elif order[char1] > order[char2]:
            return True
    return False

def sliding_window_optimal(string, k):
    window = deque(i for i in string[:k])
    min_substring = window.copy()

    for i in range(k, len(string)):
        window.append(string[i])
        window.popleft()
        if  judge(min_substring, window):
            min_substring = window.copy()
    return ''.join(min_substring)


def get_fastq(fq_file):
    my_dir=defaultdict(list)
    with open(fq_file) as file1:
        for i,line in enumerate(file1):
            if i % 4==0:
                name=line.split()[0]
            if i % 4==1:
                read=line.strip().replace("C", "T")
                my_dir[read].append(name)

    result_dir={}
    dup_dir={}
    for key,value in my_dir.items():
        if len(value)!=1:
            dup_dir[value[0]]=value[1:]
            result_dir[key]=value[0]
        else:
            result_dir[key] = value[0]
    return result_dir,dup_dir




def get_kmer_dir(fq_dir,start=0,w=30 ,k=26):
    my_dir = defaultdict(list)
    for key,value in fq_dir.items():
        my_dir[sliding_window_optimal(key[start:start + w], k)].append(key)
    return my_dir



def get_pair_with_barcode(mc_dir,c_dir,mc_fq,c_fq,dis):
    result_list = []
    for key in mc_dir.keys():
        if key in c_dir:
            list1 = mc_dir[key]
            list2 = c_dir[key]
            if len(list1) == len(list2) == 1:
                if Levenshtein.distance(list1[0], list2[0]) < dis:
                    if is_substring_with_one_error(list1[0], list2[0]):
                        result_list.append((mc_fq[list1[0]], c_fq[list2[0]]))
            else:
                temp_list = list1 + list2
                len1 = len(list1)
                temp_dir = find_dup_base_distance(temp_list,dis)
                if len(temp_dir) != 0:
                    for temp_key in temp_dir.keys():
                        temp_read1 = temp_list[int(temp_key)]
                        for temp_index in temp_dir[temp_key]:
                            temp_read2 = temp_list[temp_index]
                            if int(temp_key) < len1 and temp_index >= len1:
                                if is_substring_with_one_error(temp_read1, temp_read2):
                                    result_list.append((mc_fq[temp_read1], c_fq[temp_read2]))
    return  result_list







def dup_base_kmer(mC_file,c_file,w_file,pos,dis=4):
    mc_fq, mc1 = get_fastq(mC_file)
    mc_dir = get_kmer_dir(mc_fq, pos)

    c_fq, c1 = get_fastq(c_file)
    c_dir = get_kmer_dir(c_fq, pos)

    my_list = get_pair_with_barcode(mc_dir, c_dir, mc_fq, c_fq, dis)
    with open(w_file,"w")as file:
        for i in my_list:
            print(i[0]+" "+i[1],file=file)








if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")

    parser.add_argument("mC_file", type=str, help="Input mC file")
    parser.add_argument("c_file", type=str, help="Input c file")
    parser.add_argument("w_file", type=str, help="Output w file")
    parser.add_argument("pos", type=int, help="Position parameter")
    parser.add_argument("--dis", type=int, default=4, help="Distance parameter (default: 4)")

    args = parser.parse_args()

    mC_file = args.mC_file
    c_file = args.c_file
    w_file = args.w_file
    pos = args.pos
    dis = args.dis

    # 调用函数并传入参数
    dup_base_kmer(mC_file, c_file, w_file, pos, dis)










