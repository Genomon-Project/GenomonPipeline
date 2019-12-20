#! /usr/bin/env python

import sys, re

input_file = sys.argv[1]

re_chr_dir = re.compile('\d+\. ([^ \t\n\r\f\v,]+)-([^ \t\n\r\f\v,]+) ([fr])([fr])')

temp_chr1 = ""
temp_chr2 = ""
temp_dir1 = ""
temp_dir2 = ""
temp_gene1 = ""
temp_gene2 = ""
temp_pos1 = ""
temp_pos2 = ""
temp_read = ""
temp_pair = ""
read_count = 0
print_on = 0

with open(input_file, 'r') as hin:
    for line in hin:
        read_count = read_count + 1
        line = line.rstrip('\n')
        re_group = re_chr_dir.match(line)
        if re_group is not None:
            temp_chr1, temp_chr2, temp_dir1_raw, temp_dir2_raw = re_group.group(1), re_group.group(2), re_group.group(3), re_group.group(4)
            temp_dir1 = '+' if temp_dir1_raw == 'f' else '-'
            temp_dir2 = '+' if temp_dir2_raw == 'r' else '-'

        line = line.replace('<TR>', '')
        line = line.replace('<TD ALIGN="LEFT">', '')
        line = line.replace('<TD ALIGN="RIGHT">', '')
        line = line.replace('</TD>', '')
        line = re.sub(r'<a href=\"\#read_\d+\">', '', line)
        line = re.sub(r'<a href=\"\#pair_\d+\">', '', line)
        line = line.replace('</a>', '')

        if line.startswith('<a href="#fusion'):
            read_count = 0
            print_on = 1
            continue
    
        if read_count == 1:
            temp_gene1 = line
            continue

        if read_count == 3: 
            temp_pos1 = line
            continue

        if read_count == 4:  
            temp_gene2 = line
            continue

        if read_count == 6:  
            temp_pos2 = line
            continue

        if read_count == 7:
            temp_read = line
            continue

        if read_count == 8:
            temp_pair = line
            if print_on == 1:
                print '\t'.join([temp_chr1, temp_pos1, temp_dir1, temp_chr2, temp_pos2, temp_dir2, temp_read, temp_pair, temp_gene1, temp_gene2])
                print_on = 0

"""        
<P><P><P><BR>
18. X-X ff
<TABLE CELLPADDING=3 BORDER="1">
<TR><TD ALIGN="LEFT"><a href="#fusion_21">MCF7</a></TD>
<TD ALIGN="LEFT">TXLNG</TD>
<TD ALIGN="LEFT">X</TD>
<TD ALIGN="RIGHT">16804711</TD>
<TD ALIGN="LEFT">SYAP1</TD>
<TD ALIGN="LEFT">X</TD>
<TD ALIGN="RIGHT">16753349</TD>
<TD ALIGN="RIGHT"><a href="#read_21">20</a></TD>
<TD ALIGN="RIGHT"><a href="#pair_21">2</a></TD>
<TD ALIGN="RIGHT">16</TD>
</TR>

<H1><A NAME="detail"></A><BR>table description</H1>
"""

