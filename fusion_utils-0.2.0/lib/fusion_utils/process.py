#! /usr/bin/env python

import re

ReReg = re.compile(r'([^ \t\n\r\f\v,]+):([\+\-])(\d+)\-([^ \t\n\r\f\v,]+):([\+\-])(\d+)')

def convert_to_bedpe(input_file, output_file, margin_major, margin_minor, method):

    if method not in ["fusionfusion", "fusionfusion_part", "genomonSV", "star_fusion", "genomon_fusion", "mapsplice2", "tophat_fusion"]:
        raise ValueError("the argument method should be fusionfusion, genomonSV, star_fusion, genomon_fusion, mapsplice2, tophat_fusion")

    hIN = open(input_file, 'r')
    hOUT = open(output_file, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        if F[0] == "#fusion_name" and method == "star_fusion": continue

        chr1, pos1, dir1, chr2, pos2, dir2 = get_position(F, method)
        start1, end1, start2, end2 = pos1, pos1, pos2, pos2
        ID = chr1 + ':' + dir1 + pos1 + '-' + chr2 + ':' + dir2 + pos2

        chr1 = chr1.replace("chr", "")
        chr2 = chr2.replace("chr", "")
       
        read_num = get_read_num(F, method)

        if dir1 == '+':
            start1 = str(int(start1) - int(margin_minor))
            end1 = str(int(end1) + int(margin_major))
        else:
            start1 = str(int(start1) - int(margin_major))
            end1 = str(int(end1) + int(margin_minor))

        if dir2 == '+':
            start2 = str(int(start2) - int(margin_minor))
            end2 = str(int(end2) + int(margin_major))
        else:
            start2 = str(int(start2) - int(margin_major))
            end2 = str(int(end2) + int(margin_minor))

        print >> hOUT, '\t'.join([chr1, start1, end1, chr2, start2, end2, ID, str(read_num), dir1, dir2])

    hIN.close()
    hOUT.close()


def get_position(F, method):

    if method in ["fusionfusion", "fusionfusion_part", "genomonSV", "tophat_fusion"]:
        chr1, pos1, dir1, chr2, pos2, dir2 = F[0], F[1], F[2], F[3], F[4], F[5]
    elif method == "star_fusion":
        chr1, pos1, dir1 = F[4].split(':')
        chr2, pos2, dir2_flip = F[7].split(':')
        dir2 = '+' if dir2_flip == '-' else '-'
    elif method == "genomon_fusion":
        keyMatch = ReReg.match(F[0])
        chr1, dir1, pos1, chr2, dir2, pos2 = keyMatch.group(1), keyMatch.group(2), keyMatch.group(3), keyMatch.group(4), keyMatch.group(5), keyMatch.group(6)
    elif method == "mapsplice2":
        chr1, chr2 = F[0].split('~')
        pos1, pos2 = F[1], F[2]
        dir1, dir2 = F[5][0], '+' if F[5][1] == '-' else '-'
    
    if chr1 > chr2 or chr1 == chr2 and int(pos1) > int(pos2):
        chr1, chr2, pos1, pos2, dir1, dir2 = chr2, chr1, pos2, pos1, dir2, dir1

    return [chr1, pos1, dir1, chr2, pos2, dir2]


def get_read_num(F, method):

    read_num = 0
    if method == "fusionfusion":
        read_num = 0
        for i in range(11, len(F)):
            if F[i] != "---" and int(F[i]) > read_num:
                read_num = int(F[i])

    elif method == "fusionfusion_part":
        read_num = F[7]
    elif method == "genomonSV":
        read_num = F[13]
    elif method == "star_fusion":
        read_num = F[1]
    elif method == "genomon_fusion":
        read_num = F[20]
    elif method == "mapsplice2":
        read_num = F[21]
    elif method == "tophat_fusion":
        read_num = F[6]

    return int(read_num)
