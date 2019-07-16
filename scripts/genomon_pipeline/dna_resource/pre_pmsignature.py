#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Res_PrePmsignature(st.Stage_task):

    task_name = "pre_pmsignature"

    script_template = """
#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv
set -eu
set -o pipefail

if [ -e {output_file}.tmp ]; then
    rm {output_file}.tmp
fi

for input_file in {input_files}
do
    if [ `cat $input_file | grep -P '^[^#]' | wc -l` -gt 1 ]; then
        cut -s -f 1,2,3,5,6 $input_file | tail -n +2 | grep -P '^[^\t]+\t([Cc]hr)?[0-9XY]' > {output_file}.cut
        cut -s -f 1 {output_file}.cut > {output_file}.cut1
        cut -s -f 2 {output_file}.cut | sed "s/^/chr/" | sed -e "s/^chr[Cc]hr/chr/g" > {output_file}.cut2
        cut -s -f 3,4,5 {output_file}.cut > {output_file}.cut3
        paste {output_file}.cut1 {output_file}.cut2 {output_file}.cut3 >> {output_file}.tmp
        rm {output_file}.cut {output_file}.cut1 {output_file}.cut2 {output_file}.cut3
    fi
done

if [ -e {output_file}.tmp ]; then
    mv {output_file}.tmp {output_file}
else
    touch {output_file}
fi
"""
 
    def __init__(self, qsub_option, script_dir):
        super(Res_PrePmsignature, self).__init__(qsub_option, script_dir)

