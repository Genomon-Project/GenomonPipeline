#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class SV_filt(st.Stage_task):

    task_name = "sv_filt"

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

GenomonSV filt {input_bam} {output_prefix} {reference_genome} {param} || exit $?

mv {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.txt.tmp || exit $?

echo -e "{meta_info}" > {output_prefix}.genomonSV.result.txt || exit $?

cat {output_prefix}.genomonSV.result.txt.tmp >> {output_prefix}.genomonSV.result.txt || exit $?

rm -rf {output_prefix}.genomonSV.result.txt.tmp
 
sv_utils filter {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.filt.txt.tmp {sv_utils_param} || exit $?

mv {output_prefix}.genomonSV.result.filt.txt.tmp {output_prefix}.genomonSV.result.filt.txt
"""

    def __init__(self, qsub_option, script_dir):
        super(SV_filt, self).__init__(qsub_option, script_dir)


