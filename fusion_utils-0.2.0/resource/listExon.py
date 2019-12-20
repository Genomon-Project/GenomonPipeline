#! /usr/bin/env python

import sys, gzip

inputFile = sys.argv[1]
gene_type = sys.argv[2] # ref: refGene.txt.gz, ens: ensGene.txt.gz


with gzip.open(inputFile, 'r') as hin:
    for line in hin:

        F = line.rstrip('\n').split('\t')

        chr = F[2]
        gene_id = F[1]
        exon_starts = F[9].split(',')
        exon_ends = F[10].split(',')
        strand = F[3]
        exonNum = int(F[8])
        gene = F[1]
        symbol = F[12]

        size = 0
        for i in range(len(exon_starts) - 1):
            size = size + int(exon_ends[i]) - int(exon_starts[i])

        gene_print_name = "---"
        if gene_type == "ref":
            gene_print_name = symbol
        elif gene_type == "ens":
            gene_print_name = gene_id
        else:
            print >> sys.stderr, "The 2nd argument should be ref or ens"
            sys.exit(1)
 
        for i in range(0, len(exon_starts) - 1):
            key = chr + '\t' + exon_starts[i] + '\t' + exon_ends[i]
            if strand == "+":
                print key + '\t' + gene_print_name + '\t' + str(size) + '\t' + "+"
            else:
                print key + '\t' + gene_print_name + '\t' + str(size) + '\t' + "-"



