#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Fusion_count(st.Stage_task):

    task_name = "fusion_count"

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
export LD_LIBRARY_PATH=/usr/local/lib

chimera_utils count {additional_params} {chimeric_sam} {output}
"""

    def __init__(self, qsub_option, script_dir):
        super(Fusion_count, self).__init__(qsub_option, script_dir)
