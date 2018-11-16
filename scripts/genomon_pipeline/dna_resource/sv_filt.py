#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class SV_filt(Stage_task):

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

# set python environment
export PYTHONHOME={pythonhome}
blat_home={blat}
export PATH=${{blat_home%/*}}:{htslib}:$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

{genomon_sv} filt {input_bam} {output_prefix} {reference_genome} {param} || exit $?

mv {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.txt.tmp || exit $?

echo -e "{meta_info}" > {output_prefix}.genomonSV.result.txt || exit $?

cat {output_prefix}.genomonSV.result.txt.tmp >> {output_prefix}.genomonSV.result.txt || exit $?

rm -rf {output_prefix}.genomonSV.result.txt.tmp
 
{sv_utils} filter {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.filt.txt.tmp {sv_utils_param} || exit $?

mv {output_prefix}.genomonSV.result.filt.txt.tmp {output_prefix}.genomonSV.result.filt.txt
"""

    def __init__(self, qsub_option, script_dir):
        super(SV_filt, self).__init__(qsub_option, script_dir)


