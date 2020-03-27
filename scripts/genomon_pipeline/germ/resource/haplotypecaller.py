#! /usr/bin/env python

import os
import genomon_pipeline.core.stage_task_abc as stage_task

class Haplotypecaller(stage_task.Stage_task):
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

/usr/bin/java \\
  -XX:-UseContainerSupport \\
  -Xmx30g \\
  -jar /tools/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar HaplotypeCaller \\
  -I={INPUT_CRAM} \\
  -O={OUTPUT_VCF} \\
  -R={REFERENCE} \\
  --native-pair-hmm-threads=$(nproc) \\
  --TMP_DIR=$(dirname {OUTPUT_VCF})
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(input_bams, genomon_conf, run_conf, sample_conf):
    
    STAGE_NAME = "gatk-haplotypecaller-parabrics-compatible"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.get(CONF_SECTION, "image"),
        "qsub_option": genomon_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": genomon_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Haplotypecaller(params)
    
    output_files = []
    for sample in sample_conf.haplotype_call:
        output_vcf = "haplotypecaller/%s/%s.gatk-hc.vcf" % (sample, sample)
        output_files.append(output_vcf)
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_VCF":  "%s/%s" % (run_conf.project_root, output_vcf),
            "REFERENCE": genomon_conf.path_get(CONF_SECTION, "reference"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(genomon_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
    return output_files
