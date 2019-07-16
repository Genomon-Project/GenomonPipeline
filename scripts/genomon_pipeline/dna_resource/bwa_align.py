#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Bwa_align(st.Stage_task):

    task_name = "bwa_align"

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

tmp_num=`expr ${{SGE_TASK_ID}} - 1`
num=`printf "%04d" ${{tmp_num}}`

{bwa} mem {bwa_params} {ref_fa} {input_dir}/1_${{num}}.fastq_split {input_dir}/2_${{num}}.fastq_split > {output_dir}/{sample_name}_${{num}}.bwa.sam || exit $?

{biobambam}/bamsort index=1 level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindentonly=1 calmdnmreference={ref_fa} tmpfile={output_dir}/{sample_name}_${{num}}.sorted.bam.tmp inputformat=sam indexfilename={output_dir}/{sample_name}_${{num}}.sorted.bam.bai I={output_dir}/{sample_name}_${{num}}.bwa.sam O={output_dir}/{sample_name}_${{num}}.sorted.bam

"""

    def __init__(self, qsub_option, script_dir):
        super(Bwa_align, self).__init__(qsub_option, script_dir)


