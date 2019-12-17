#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Res_QC_Bamstats(st.Stage_task):

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
export bamstats=/tools/ICGC/bin/bam_stats.pl
export ld_library_path=/usr/local/bin
export perl5lib=/tools/ICGC/lib/perl5:/tools/ICGC/lib/perl5/x86_64-linux-gnu-thread-multi

genomon_qc bamstats {input_file} {output_file} --perl5lib $perl5lib --bamstats $bamstats
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Bamstats, self).__init__(qsub_option, script_dir)
