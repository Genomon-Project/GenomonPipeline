#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class SV_merge(Stage_task):

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

# set python environment
export PYTHONHOME={pythonhome}
export PATH={htslib}:$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

{genomon_sv} merge {control_info} {merge_output_file} {param} || exit $?  

"""

    def __init__(self, qsub_option, script_dir):
        super(SV_merge, self).__init__(qsub_option, script_dir)


