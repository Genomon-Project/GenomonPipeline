#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Genomon_expression(st.Stage_task):

    task_name = "genomon_expression"

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

genomon_expression {additional_params} {input_bam} {output_prefix}
mv {output_prefix}.sym2fpkm.txt {output_prefix}.genomonExpression.result.txt
"""

    def __init__(self, qsub_option, script_dir):
        super(Genomon_expression, self).__init__(qsub_option, script_dir)
