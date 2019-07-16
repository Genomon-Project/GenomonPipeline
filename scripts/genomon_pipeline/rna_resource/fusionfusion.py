#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Fusionfusion(st.Stage_task):

    task_name = "fusionfusion"

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

{fusionfusion} --star {chimeric_sam} --out {output_prefix} --reference_genome {ref_fa} {additional_params} || exit $?

mv {output_prefix}/star.fusion.result.txt {output_prefix}/{sample}.star.fusion.result.txt || exit $?
mv {output_prefix}/fusion_fusion.result.txt {output_prefix}/{sample}.genomonFusion.result.txt || exit $?

{fusion_utils} filt {output_prefix}/{sample}.genomonFusion.result.txt {output_prefix}/{sample}.fusion.fusion.result.filt.txt {filt_params} || exit $?
mv {output_prefix}/{sample}.fusion.fusion.result.filt.txt {output_prefix}/{sample}.genomonFusion.result.filt.txt
"""

    def __init__(self, qsub_option, script_dir):
        super(Fusionfusion, self).__init__(qsub_option, script_dir)
