#! /usr/bin/env python

import os
import genomon_pipeline.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "cram/{sample}/{sample}.markdup.cram"
BAM_POSTFIX = ".markdup.cram"
BAI_POSTFIX = ".markdup.cram.crai"
    
class Fq2cram(stage_task.Stage_task):
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

work_dir=$(dirname {OUTPUT_CRAM})

/tools/bwa-0.7.15/bwa mem \\
  -t $(nproc) \\
  -K 10000000 \\
  -T 0 \\
  -R "@RG\\tID:{SAMPLE_NAME}\\tPL:na\\tLB:na\\tSM:{SAMPLE_NAME}\\tPU:{SAMPLE_NAME}" \\
  {REFERENCE} \\
  {INPUT_FASTQ_1} \\
  {INPUT_FASTQ_2} \\
| /usr/bin/java \\
  -XX:-UseContainerSupport \\
  -Xmx30g \\
  -jar /tools/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar SortSam \\
  --MAX_RECORDS_IN_RAM=5000000 \\
  -I=/dev/stdin \\
  -O=${{work_dir}}/{SAMPLE_NAME}.bam \\
  --SORT_ORDER=coordinate \\
  --TMP_DIR=${{work_dir}}

/usr/bin/java \\
  -XX:-UseContainerSupport \\
  -Xmx30g \\
  -jar /tools/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar MarkDuplicates \\
  -I=${{work_dir}}/{SAMPLE_NAME}.bam \\
  -O=${{work_dir}}/{SAMPLE_NAME}.markdup.bam \\
  -M={OUTPUT_MARKDUP_METRICS} \\
  --TMP_DIR=${{work_dir}}

rm -f ${{work_dir}}/{SAMPLE_NAME}.bam

/tools/samtools-1.9/samtools view \\
  -@ $(nproc) \\
  -C \\
  -T {REFERENCE} \\
  -o {OUTPUT_CRAM} \\
  ${{work_dir}}/{SAMPLE_NAME}.markdup.bam

/tools/samtools-1.9/samtools index \\
  -@ $(nproc) \\
  {OUTPUT_CRAM}

rm -f ${{work_dir}}/{SAMPLE_NAME}.markdup.bam
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(genomon_conf, run_conf, sample_conf):

    STAGE_NAME = "bwa-alignment-parabrics-compatible"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": genomon_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": genomon_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Fq2cram(params)
     
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/cram/%s" % (run_conf.project_root, sample)
        output_bams[sample] = "%s/%s.markdup.cram" % (output_dir, sample)
        
        if len(sample_conf.fastq[sample][0]) == 1:
            fastq1 = sample_conf.fastq[sample][0][0]
        else:
            fastq1 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        if len(sample_conf.fastq[sample][1]) == 1:
            fastq2 = sample_conf.fastq[sample][1][0]
        else:
            fastq2 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        arguments = {
            "SAMPLE_NAME": sample,
            "INPUT_FASTQ_1": fastq1,
            "INPUT_FASTQ_2": fastq2,
            "OUTPUT_CRAM": output_bams[sample],
            "OUTPUT_MARKDUP_METRICS": "%s/%s.markdup.metrics" % (output_dir, sample),
            "REFERENCE": genomon_conf.path_get(CONF_SECTION, "reference"),
        }
        
        singularity_bind = [
            run_conf.project_root,
            os.path.dirname(genomon_conf.path_get(CONF_SECTION, "reference")),
        ] + sample_conf.fastq_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams
