#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Res_Pmsignature(st.Stage_task):

    task_name = "pmsignature"

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

export PATH={r_path}:$PATH
export LD_LIBRARY_PATH={r_ld_library_path}:$LD_LIBRARY_PATH
export R_LIBS={r_libs}
export R_PATH={r_path}

sig_list=({sig_list})
sig_num=${{sig_list[$SGE_TASK_ID-1]}}

{command} 
"""

    ind_template = """
$R_PATH/R --vanilla --slave --args {inputfile} {outputdir}/ind.$sig_num.Rdata $sig_num {trdirflag} {trialnum} {bgflag} {bs_genome} {txdb_transcript} < {script_path}/pmsignature/run_pmsignature_ind.R
if [ $? -ne 0 ]
then
    echo pmsignature terminated abnormally.
    echo '{{"id":[],"ref":[],"alt":[],"strand":[],"mutation":[]}}' > {outputdir}/pmsignature.ind.result.$sig_num.json
    exit 0
fi
$R_PATH/R --vanilla --slave --args {outputdir}/ind.$sig_num.Rdata {outputdir}/pmsignature.ind.result.$sig_num.json < {script_path}/pmsignature/convert_toJson_ind.R
"""

    full_template = """
$R_PATH/R --vanilla --slave --args {inputfile} {outputdir}/full.$sig_num.Rdata $sig_num {trdirflag} {trialnum} {bgflag} {bs_genome} {txdb_transcript} < {script_path}/pmsignature/run_pmsignature_full.R
if [ $? -ne 0 ]
then
    echo pmsignature terminated abnormally.
    echo '{{"id":[],"signature":[],"mutation":[]}}' > {outputdir}/pmsignature.full.result.$sig_num.json
    exit 0
fi
$R_PATH/R --vanilla --slave --args {outputdir}/full.$sig_num.Rdata {outputdir}/pmsignature.full.result.$sig_num.json < {script_path}/pmsignature/convert_toJson_full.R
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_Pmsignature, self).__init__(qsub_option, script_dir)

