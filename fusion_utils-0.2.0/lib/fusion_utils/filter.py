#! /usr/bin/env python

def filter_star_fusion(input_file, output_file, thres):

    hIN = open(input_file, 'r')
    hOUT = open(output_file, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        if F[0] == "#fusion_name":
            print >> hOUT, '\t'.join(F)
            continue 

        if int(F[1]) < int(thres): continue

        if F[4].startswith("MT"): continue
        if F[4].startswith("GL0"): continue 
        if F[7].startswith("MT"): continue
        if F[7].startswith("GL0"): continue

        print >> hOUT, '\t'.join(F)

    hIN.close()
    hOUT.close()


def filter_genomon_fusion(input_file, output_file, thres):

    hIN = open(input_file, 'r')
    hOUT = open(output_file, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        if "hs37d5" in F[0]: continue

        non_gene_flag = 0
        # if F[7] == "---" and F[8] == "---": non_gene_flag = 1
        if ("NM_" not in F[7] and "NR_" not in F[7]) or ("NM_" not in F[8] and "NR_" not in F[8]): non_gene_flag = 1

        non_junc_flag = 0
        # if F[9] == "---" and F[10] == "---": non_junc_flag = 1
        if "NM_" not in F[9] and "NR_" not in F[9] and "NM_" not in F[10] and "NR_" not in F[10]: non_junc_flag = 1

        if non_gene_flag == 1 and non_junc_flag == 1: continue

        if int(F[20]) < thres: continue

        print >> hOUT, '\t'.join(F)

    hIN.close()
    hOUT.close()


def filter_fusionfusion(input_file, output_file, thres):

    hIN = open(input_file, 'r')
    hOUT = open(output_file, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        if F[0] == "hs37d5" or F[3] == "hs37d5": continue

        non_gene_flag = 0
        if F[8] == "---" or F[10] == "---": non_gene_flag = 1

        non_junc_flag = 0
        if F[9] == "---" and F[11] == "---": non_junc_flag = 1

        if non_gene_flag == 1 and non_junc_flag == 1: continue

        if int(F[7]) < thres: continue

        print >> hOUT, '\t'.join(F)

    hIN.close()
    hOUT.close()

