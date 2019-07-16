#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Intron_retention(st.Stage_task):

    task_name = "intron_retention"

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

# set python environment
export PYTHONHOME={pythonhome}
bedtools_home={bedtools}
export PATH=${{bedtools_home%/*}}:{htslib}:$PYTHONHOME/bin:$PATH
export PYTHONPATH={pythonpath}

{intron_retention_utils} simple_count {additional_params} {input_bam} {output_prefix}.ir_simple_count.txt
mv {output_prefix}.ir_simple_count.txt {output_prefix}.genomonIR.result.txt
"""

    def __init__(self, qsub_option, script_dir):
        super(Intron_retention, self).__init__(qsub_option, script_dir)

