#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_QC_Bamstats(Stage_task):

    task_name = "qc_bamstats"

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
export PYTHONPATH={pythonpath}
export PATH=$PYTHONHOME/bin:{perl}:$PATH
export PERL5LIB={perl5lib}

{genomon_qc} bamstats {input_file} {output_file} --bamstats {bamstats}
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Bamstats, self).__init__(qsub_option, script_dir)
