#! /usr/bin/env python

import os
import genomon_pipeline.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "bam/{sample}/{sample}.markdup.bam"
BAM_POSTFIX = ".markdup.bam"
BAI_POSTFIX = ".markdup.bam.bai"
        
class Bwa_align(stage_task.Stage_task):
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

export LD_LIBRARY_PATH=/usr/local/lib
OUTPUT_PREF={OUTPUT_DIR}/{SAMPLE}
mkdir -p {OUTPUT_DIR}

/tools/bwa-0.7.17/bwa mem \
  {BWA_OPTION} \
  {REFERENCE} \
  {FASTQ1} \
  {FASTQ2} \
| /usr/local/bin/bamsort \
  {BAMSORT_OPTION} \
  calmdnmreference={REFERENCE} \
  inputformat=sam \
  indexfilename=${{OUTPUT_PREF}}.sorted.bam.bai \
  O=${{OUTPUT_PREF}}.sorted.bam

/usr/local/bin/bammarkduplicates \
  {BAMMARKDUP_OPTION} \
  M=${{OUTPUT_PREF}}.metrics \
  I=${{OUTPUT_PREF}}.sorted.bam \
  O=${{OUTPUT_PREF}}.markdup.bam

rm ${{OUTPUT_PREF}}.sorted.bam
rm ${{OUTPUT_PREF}}.sorted.bam.bai
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(genomon_conf, run_conf, sample_conf):

    STAGE_NAME = "bwa_alignment"
    CONF_SECTION = "bwa_alignment"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": genomon_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": genomon_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Bwa_align(params)
     
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/bam/%s" % (run_conf.project_root, sample)
        output_bams[sample] = "%s/%s.markdup.bam" % (output_dir, sample)
        
        if len(sample_conf.fastq[sample][0]) == 1:
            fastq1 = sample_conf.fastq[sample][0][0]
        else:
            fastq1 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        if len(sample_conf.fastq[sample][1]) == 1:
            fastq2 = sample_conf.fastq[sample][1][0]
        else:
            fastq2 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        arguments = {
            "SAMPLE": sample,
            "INPUT_BAM": sample_conf.fastq[sample],
            "FASTQ1": fastq1,
            "FASTQ2": fastq2,
            "OUTPUT_DIR": output_dir,
            "REFERENCE": genomon_conf.path_get(CONF_SECTION, "bwa_reference"),
            "BWA_OPTION": genomon_conf.get(CONF_SECTION, "bwa_option"),
            "BAMSORT_OPTION": genomon_conf.get(CONF_SECTION, "bamsort_option"),
            "BAMMARKDUP_OPTION": genomon_conf.get(CONF_SECTION, "bammarkduplicates_option")
        }
        
        singularity_bind = [
            run_conf.project_root,
            os.path.dirname(genomon_conf.path_get(CONF_SECTION, "bwa_reference")),
        ] + sample_conf.fastq_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams
