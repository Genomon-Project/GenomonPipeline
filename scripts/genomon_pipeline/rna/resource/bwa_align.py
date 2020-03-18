#! /usr/bin/env python

import genomon_pipeline.core.stage_task_abc as stage_task

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
REFERENCE={REFERENCE_DIR}/{REFERENCE_FILE}
mkdir -p {OUTPUT_DIR}

/tools/bwa-0.7.17/bwa mem \
  {BWA_OPTION} \
  ${{REFERENCE}} \
  {FASTQ1} \
  {FASTQ2} \
| /usr/local/bin/bamsort \
  {BAMSORT_OPTION} \
  calmdnmreference=${{REFERENCE}} \
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
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": "bwa_alignment",
        "image": genomon_conf.get("bwa_alignment", "image"),
        "qsub_option": genomon_conf.get("bwa_alignment", "qsub_option"),
        "singularity_option": genomon_conf.get("bwa_alignment", "singularity_option")
    }
    stage_class = Bwa_align(params)
    
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/bam/%s" % (run_conf.project_root, sample)
        output_bams[sample] = "%s/%s.markdup.bam" % (output_dir, sample)
        
        arguments = {
            "SAMPLE": sample,
            "INPUT_BAM": sample_conf.fastq[sample],
            "FASTQ1": sample_conf.fastq[sample][0][0],
            "FASTQ2": sample_conf.fastq[sample][1][0],
            "OUTPUT_DIR": output_dir,
            "REFERENCE_DIR": genomon_conf.get("bwa_alignment", "bwa_reference_dir"),
            "REFERENCE_FILE": genomon_conf.get("bwa_alignment", "bwa_reference_file"),
            "BWA_OPTION": genomon_conf.get("bwa_alignment", "bwa_option"),
            "BAMSORT_OPTION": genomon_conf.get("bwa_alignment", "bamsort_option"),
            "BAMMARKDUP_OPTION": genomon_conf.get("bwa_alignment", "bammarkduplicates_option")
        }
        
        singularity_bind = [
            run_conf.project_root,
            genomon_conf.get("bwa_alignment", "bwa_reference_dir"),
        ] + sample_conf.fastq_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams
