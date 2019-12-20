#! /usr/bin/env python

import sys, gzip

inputFile = sys.argv[1]
gene_type = sys.argv[2] # ref: refGene.txt.gz, ens: ensGene.txt.gz

with gzip.open(inputFile, 'r') as hin:
    for line in hin:

        F = line.rstrip('\n').split('\t')

        chr = F[2]
        gene_id = F[1]
        gene_start = F[4]
        gene_end = F[5]
        strand = F[3]
        symbol = F[12]
        exon_starts = F[9].split(',')
        exon_ends = F[10].split(',')

        size = 0
        for i in range(len(exon_starts) - 1):
            size = size + int(exon_ends[i]) - int(exon_starts[i])

        key = chr + '\t' + gene_start + '\t' + gene_end
        if gene_type == "ref": 
            print key + '\t' + symbol + '\t' + str(size) + '\t' + strand
        elif gene_type == "ens":
            print key + '\t' + gene_id + '\t' + str(size) + '\t' + strand
        else:
            print >> sys.stderr, "The 2nd argument should be ref or ens"
            sys.exit(1)


