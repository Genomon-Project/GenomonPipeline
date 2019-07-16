#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Res_PA_Plot(st.Stage_task):

    task_name = "paplot"

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
set -o pipefail

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export PYTHONPATH={pythonpath}

if test "{inputs_qc}" != ""; then
    {paplot} qc "{inputs_qc}" {output_dir} {title} --config_file {config_file} --remarks "{remarks}" --title 'QC graphs' --overview 'Quality Control of bam.' --ellipsis qc
fi
if test "{inputs_sv}" != ""; then
    {paplot} ca "{inputs_sv}" {output_dir} {title} --config_file {config_file} --remarks "{remarks}" --title 'Fusion graphs' --overview 'Fusion.' --ellipsis fusion
fi
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PA_Plot, self).__init__(qsub_option, script_dir)

        
