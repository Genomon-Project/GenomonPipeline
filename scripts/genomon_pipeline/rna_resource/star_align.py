#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Star_align(st.Stage_task):

    task_name = "star_align"

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

export LD_LIBRARY_PATH=/usr/local/lib

/usr/local/bin/STAR --genomeDir {star_genome} --readFilesIn {fastq1} {fastq2} --outFileNamePrefix {out_prefix} {additional_params} 

/usr/local/bin/samtools sort -T {out_prefix}Aligned.sortedByCoord.out {samtools_sort_params} {out_prefix}Aligned.out.bam -O bam > {out_prefix}Aligned.sortedByCoord.out.bam 

/usr/local/bin/samtools index {out_prefix}Aligned.sortedByCoord.out.bam 

rm -rf {out_prefix}Aligned.out.bam
"""

    def __init__(self, qsub_option, script_dir):
        super(Star_align, self).__init__(qsub_option, script_dir)
