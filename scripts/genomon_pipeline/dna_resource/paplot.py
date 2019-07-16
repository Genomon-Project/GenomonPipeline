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

{command} 
"""

    qc_template = """
{paplot} qc '{inputs}' {output_dir} {title} --config_file {config_file} --title 'QC graphs' --overview 'Quality Control of bam.' --ellipsis qc
"""
    sv_template = """
{paplot} ca '{inputs}' {output_dir} {title} --config_file {config_file} --title 'SV graphs' --overview 'Structural Variation.' --ellipsis sv 
"""
    mutation_template = """
if test '{annovar}' == 'True'; then
    {paplot} mutation '{inputs}' {output_dir} {title} --config_file {config_file} --title 'Mutation matrix' --overview 'Gene-sample mutational profiles.' --ellipsis mutation
else
    echo 'paplot: [annotation] active_annovar_flag = False in genomon_conf_file, skip mutation-matrix.'
fi
"""
    full_template = """
{paplot} signature {input} {output_dir} {title} --config_file {config_file} --title 'Mutational Signature' --overview 'Pmsignature type=full.' --ellipsis full
"""
    ind_template = """
{paplot} pmsignature {input} {output_dir} {title} --config_file {config_file} --title 'pmsignature' --overview 'Pmsignature type=ind.' --ellipsis ind
"""

    index_template = """
{paplot} index {output_dir} --config_file {config_file} --remarks '{remarks}'
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PA_Plot, self).__init__(qsub_option, script_dir)

        
