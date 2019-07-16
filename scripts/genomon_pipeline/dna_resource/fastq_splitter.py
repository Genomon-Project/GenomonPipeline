#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Fastq_splitter(st.Stage_task):

    task_name = "fastq_splitter"

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
set -o pipefail

to_val=`ls {target_dir}/*_${{SGE_TASK_ID}}{ext} | wc -l`
input_files=""
for i in `seq 1 ${{to_val}}`; do
    input_files="${{input_files}} {target_dir}/${{i}}_${{SGE_TASK_ID}}{ext}"
done    

if [ "_{fastq_filter}" = "_True" ]; then

    if [ "_{ext}" = "_.gz" ]; then
        gzip -dc $input_files | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    else
        cat $input_files | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    fi
else
    if [ "_{ext}" = "_.gz" ]; then
        gzip -dc $input_files | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    else
        cat $input_files | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    fi
fi

ls -1 {target_dir}/${{SGE_TASK_ID}}_[0-9][0-9][0-9][0-9] | while read filename; do
    mv $filename $filename.fastq_split || exit $?
done


"""

    def __init__(self, qsub_option, script_dir):
        super(Fastq_splitter, self).__init__(qsub_option, script_dir)


