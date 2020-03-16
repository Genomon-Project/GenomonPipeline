#! /usr/bin/env python

import genomon_pipeline.stage_task_abc as stage_task

class Sv_dummy(stage_task.Stage_task):
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

mkdir -p {OUTPUT_DIR}
samtools view {INPUT_BAM} | wc -l > {OUTPUT_DIR}/{SAMPLE}.txt
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(output_bams, genomon_conf, run_conf, sample_conf):
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": "sv_dummy",
        "image": genomon_conf.get("mutation_dummy", "image"),
        "qsub_option": genomon_conf.get("mutation_dummy", "qsub_option"),
        "singularity_option": genomon_conf.get("mutation_dummy", "singularity_option")
    }
    stage_class = Sv_dummy(params)
    
    for sample in output_bams:
        output_dir = "%s/sv/%s" % (run_conf.project_root, sample)
    
        arguments = {
            "SAMPLE": sample,
            "INPUT_BAM": output_bams[sample],
            "OUTPUT_DIR": output_dir,
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]

        stage_class.write_script(arguments, singularity_bind, sample = sample)
