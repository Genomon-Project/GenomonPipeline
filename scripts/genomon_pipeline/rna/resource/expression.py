#! /usr/bin/env python

import genomon_pipeline.core.stage_task_abc as stage_task

class Expression(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """
#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -o errexit
set -o nounset
set -o pipefail
set -x

OUTPUT_PREF={OUTPUT_DIR}/{SAMPLE}
mkdir -p {OUTPUT_DIR}

featureCounts -T 4 -p -a {GTF} -O -B -C -o ${{OUTPUT_PREF}}.txt {INPUT_BAM}
python /tools/simple_exp/proc_fc.py ${{OUTPUT_PREF}}.txt ${{OUTPUT_PREF}}.txt.summary {GTF} > ${{OUTPUT_PREF}}.txt.fpkm
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(input_bams, genomon_conf, run_conf, sample_conf):
    STAGE_NAME = "expression"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": genomon_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": genomon_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Expression(params)
    
    output_files = []
    for sample in sample_conf.expression:
        output_dir = "%s/expression/%s" % (run_conf.project_root, sample)
        output_files.append("expression/{sample}/{sample}.txt.fpkm".format(sample = sample))
        
        arguments = {
            "SAMPLE": sample,
            "INPUT_BAM": input_bams[sample],
            "OUTPUT_DIR": output_dir,
            "GTF": genomon_conf.path_get(SECTION_NAME, "gtf")
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
