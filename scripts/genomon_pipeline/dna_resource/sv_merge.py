#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class SV_merge(st.Stage_task):

    task_name = "sv_merge"

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

GenomonSV merge {control_info} {merge_output_file} {param} || exit $?  

"""

    def __init__(self, qsub_option, script_dir):
        super(SV_merge, self).__init__(qsub_option, script_dir)


