#! /usr/bin/env python

import genomon_pipeline.core.stage_task_abc as stage_task

class Star_align(stage_task.Stage_task):
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

# STAR
/usr/local/bin/STAR \
  --genomeDir {STAR_REFERENCE} \
  --readFilesIn {FASTQ1} {FASTQ2} \
  --outFileNamePrefix ${{OUTPUT_PREF}}. \
  {STAR_OPTION}

# sort
/usr/local/bin/samtools sort \
  -T ${{OUTPUT_PREF}}.Aligned.sortedByCoord.out \
  -@ 6 -m 3G \
  ${{OUTPUT_PREF}}.Aligned.out.bam \
  -O bam > ${{OUTPUT_PREF}}.Aligned.sortedByCoord.out.bam

rm ${{OUTPUT_PREF}}.Aligned.out.bam

# index
/usr/local/bin/samtools index \
  ${{OUTPUT_PREF}}.Aligned.sortedByCoord.out.bam

rm {RM_FASTQ}
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(genomon_conf, run_conf, sample_conf):
    STAGE_NAME = "star_alignment"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": genomon_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": genomon_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Star_align(params)
    
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/star/%s" % (run_conf.project_root, sample)
        output_bams[sample] = "%s/%s.Aligned.sortedByCoord.out.bam" % (output_dir, sample)

        rm_fastq = " ".join(sample_conf.fastq[sample][0])
        if len(sample_conf.fastq[sample][0]) == 1:
            fastq1 = sample_conf.fastq[sample][0][0]
        else:
            fastq1 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))
        
        fastq2 = ""
        if len(sample_conf.fastq[sample]) == 2:
            rm_fastq += " " + " ".join(sample_conf.fastq[sample][1])
            if len(sample_conf.fastq[sample][1]) == 1:
                fastq2 = sample_conf.fastq[sample][1][0]
            else:
                fastq2 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        arguments = {
            "SAMPLE": sample,
            "FASTQ1": fastq1,
            "FASTQ2": fastq2,
            "RM_FASTQ": rm_fastq,
            "OUTPUT_DIR": output_dir,
            "STAR_REFERENCE": genomon_conf.path_get(SECTION_NAME, "star_genome"),
            "STAR_OPTION": genomon_conf.get(SECTION_NAME, "star_option"),
        }
        
        singularity_bind = [
            run_conf.project_root,
            genomon_conf.get(SECTION_NAME, "star_genome"),
        ] + sample_conf.fastq_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams
