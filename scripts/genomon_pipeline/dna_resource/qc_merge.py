#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_QC_Merge(Stage_task):
    
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

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export PYTHONPATH={pythonpath}

{genomon_qc} merge {coverage_file} {bamstats_file} {output_file} --meta {meta}
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Merge, self).__init__(qsub_option, script_dir)
