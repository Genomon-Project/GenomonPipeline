#! /usr/bin/env python

import pysam

junction_margin = 5

def get_gene_info(chr, pos, ref_gene_tb, ens_gene_tb):

    # check gene annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_gene_tb.fetch(chr, int(pos) - 1, int(pos) + 1)
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1
        
    gene = [];
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            gene.append(record[3])

    # if any gene cannot be found in refGene, then search ensGene
    if len(gene) == 0:
        tabixErrorFlag = 0
        try:
            records = ens_gene_tb.fetch(chr, int(pos) - 1, int(pos) + 1)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1
            
        # for ensGene, just the longest gene is shown
        temp_length = 0
        temp_gene = ""
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if int(record[4]) > temp_length:
                    temp_gene = record[3]

            if temp_gene != "": gene.append(temp_gene)
            
    # if len(gene) == 0: gene.append("---")

    return list(set(gene))


def get_junc_info(chr, pos, ref_exon_tb, ens_exon_tb, junction_margin):

    # check exon annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_exon_tb.fetch(chr, int(pos) - junction_margin, int(pos) + junction_margin)
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1
        
    junction = []
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if abs(int(pos) - int(record[1])) < junction_margin:
                if record[5] == "+": junction.append(record[3] + ".start")
                if record[5] == "-": junction.append(record[3] + ".end")
            if abs(int(pos) - int(record[2])) < junction_margin:
                if record[5] == "+": junction.append(record[3] + ".end")
                if record[5] == "-": junction.append(record[3] + ".start")

    # if any exon-intron junction cannot be found in refGene, then search ensGene
    if len(junction) == 0: 
        tabixErrorFlag = 0
        try:
            records = ens_exon_tb.fetch(chr, int(pos) - junction_margin, int(pos) + junction_margin)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorFlag = 1
             
        # for ensGene, just the longest gene is shown
        temp_length = 0
        temp_junc = ""
        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if int(record[4]) > temp_length:
                    if abs(int(pos) - int(record[1])) < junction_margin:
                        if record[5] == "+": temp_junc = record[3] + ".start"
                        if record[5] == "-": temp_junc = record[3] + ".end"
                    if abs(int(pos) - int(record[2])) < junction_margin: 
                        if record[5] == "+": temp_junc = record[3] + ".end"
                        if record[5] == "-": temp_junc = record[3] + ".start"
                
            if temp_junc != "": junction.append(temp_junc)

                
    # if len(junction) == 0: junction.append("---")
    
    return list(set(junction))


