#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Res_QC_Merge(st.Stage_task):
    
    task_name = "qc_merge"

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

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export PYTHONPATH={pythonpath}

if [ -f {fastq_line_num_file} ]; then

    total_reads=`awk 'NR==2 {{print $15}}' {bamstats_file}`
    fastq_reads_tmp=`cat {fastq_line_num_file}`
    fastq_reads=`expr $fastq_reads_tmp / 2`

    if [ $total_reads -ne $fastq_reads ]; then
        echo "Total read count is not good for this data. BAM file: ${{total_reads}} reads. FASTQ file: ${{fastq_reads}} reads." >&2
        exit 1
    fi
fi

{genomon_qc} merge {coverage_file} {bamstats_file} {output_file} --meta "{meta}"
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Merge, self).__init__(qsub_option, script_dir)
