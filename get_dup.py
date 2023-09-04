from Bio import pairwise2
from collections import defaultdict
import Levenshtein

def is_substring_with_one_error(main_string, substring):
    start = main_string.find(substring)
    if start != -1:
        return start, start + len(substring)

    # max_allowed_errors = 1
    # alignments = pairwise2.align.localms(main_string, substring, 1, -1, -1, -1)
    #
    # for alignment in alignments:
    #     alignment_score = alignment.score
    #     if alignment_score >= len(substring) - max_allowed_errors:
    #         return alignment.start, alignment.end
    return -1

def write_file_bar(name,seq,quality,cut_pos,barcode,file_name):
    print(name + " " + barcode, file=file_name)
    print(seq[cut_pos:], file=file_name)
    print("+",file=file_name)
    print(quality[cut_pos:], file=file_name)


bar_dir_mc1=defaultdict(list)
bar_dir_c1 = defaultdict(list)

with open("test_1.fq") as file, open("test_2.fq") as file2, open("mC_bar.f1.fq", "w") as mc_bar1, open("C_bar.f1.fq",
        "w") as c_bar1, open("mC_bar.f2.fq", "w") as mc_bar2,open("C_bar.f2.fq", "w") as c_bar2, open("mC_nobar.f1.fq",
            "w") as mc_nobar1,open("mC_nobar.f2.fq", "w") as mc_nobar2 ,open("C_nobar.f1.fq",
            "w") as c_nobar1,open("C_nobar.f2.fq", "w") as c_nobar2:
    for i, line in enumerate(zip(file, file2)):
        if i % 4 == 0:
            name1 = line[0].split()[0]
            name2 = line[1].split()[0]
        if i % 4 == 1:
            seq1 = line[0].strip()
            seq2 = line[1].strip()
        if i % 4 == 3:
            quality1 = line[0].strip()
            quality2 = line[1].strip()

            add_status=""

            barcode1 = "N"
            barcode2 = "N"

            cut_pos1 = 0
            cut_pos2 = 0

            start_end1 = is_substring_with_one_error(seq2[:35], "GTGTATAAGAGA")
            start_end2 = is_substring_with_one_error(seq2[:35], "ATATATAAAAAA")

            start_end3 = is_substring_with_one_error(seq1[:35], "ATGTGTATAAGA")

            if start_end1 != -1 and start_end2 == -1:
                add_status = "mC"
                if start_end1[0] > 11:
                    barcode2 = seq2[start_end1[0]  - 12 :start_end1[0] - 4]
                cut_pos2 = start_end1[1] + 3

            elif start_end2 != -1 and start_end1 == -1:
                add_status = "C"
                if start_end2[0] > 11:
                    barcode2 = seq2[start_end2[0] - 12:start_end2[0]-4]
                cut_pos2 = start_end2[1] + 3

            else:
                if seq1[:int(len(seq1)/2)].count("C") + seq2[:int(len(seq2)/2)].count("G") > int(len(seq1)/2)+int(len(seq2)/2)*0.03:
                    add_status = "mC"
                else:
                    add_status = "C"

                pass


            if start_end3 != -1:
                if start_end3[0] > 9:
                    barcode1 = seq1[start_end3[0] - 10:start_end3[0]-2]
                cut_pos1 = start_end3[1] + 5
            else:
                pass

            if add_status=="mC" and  cut_pos1!=0 and cut_pos2!=0 :
                write_file_bar(name1,seq1,quality1,cut_pos1,barcode1,mc_bar1)
                write_file_bar(name2, seq2, quality2, cut_pos2, barcode2, mc_bar2)

            elif add_status=="C" and  cut_pos1!=0 and cut_pos2!=0:
                write_file_bar(name1, seq1, quality1, cut_pos1, barcode1, c_bar1)
                write_file_bar(name2, seq2, quality2, cut_pos2, barcode2, c_bar2)

            elif add_status=="mC":
                write_file_bar(name1, seq1, quality1, cut_pos1, barcode1, mc_nobar1)
                write_file_bar(name2, seq2, quality2, cut_pos2, barcode2, mc_nobar2)

            elif add_status=="C":
                write_file_bar(name1, seq1, quality1, cut_pos1, barcode1, c_nobar1)
                write_file_bar(name2, seq2, quality2, cut_pos2, barcode2, c_nobar2)



            if add_status=="mC" and len(barcode1) + len(barcode2)!=2:
                bar_dir_mc1[barcode1.replace("C","T")+barcode2.replace("G","A")].append(name1+" "+seq1[cut_pos1:])
            elif add_status=="C" and len(barcode1) + len(barcode2)!=2:
                bar_dir_c1[barcode1.replace("C","T")+barcode2.replace("G","A")].append(name1+" "+seq1[cut_pos1:])



yy=0
pp=0
result=[]
for key in bar_dir_mc1.keys():
    pp = pp + 1
    # if key == "TATTTTTTTACAAAAA" or key=="TGTGTGTTAAAAAATC":
    #     continue
    if key in bar_dir_c1:
        yy=yy+1
        mc_list =  bar_dir_c1[key]
        c_list =   bar_dir_mc1[key]
        for mc_seq in mc_list:
            mc_read=mc_seq.split()[1]
            for c_seq in c_list:
                c_read=c_seq.split()[1]
                len_seq=min(len(mc_read),len(c_read))
                if Levenshtein.hamming(mc_read[:len_seq].replace("C","T"),c_read[:len_seq].replace("C","T"))<4:
                    result.append((mc_seq.split()[0],c_seq.split()[0]))






    if pp %10000==0:
        print(pp)
        print(yy)

print(len(result))









