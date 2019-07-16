#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Fusion_merge(st.Stage_task):

    task_name = "fusion_merge"

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
export PATH={htslib}:$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

{chimera_utils} merge_control {additional_params} {count_list} {output}
"""

    def __init__(self, qsub_option, script_dir):
        super(Fusion_merge, self).__init__(qsub_option, script_dir)
