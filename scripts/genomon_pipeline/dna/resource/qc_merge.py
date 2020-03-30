#! /usr/bin/env python

import genomon_pipeline.core.stage_task_abc as stage_task

class Qc_merge(stage_task.Stage_task):
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
set -xv
set -eu

# set python environment
if [ -f {fastq_line_num_file} ]; then

    total_reads=`awk 'NR==2 {{print $15}}' {bamstats_file}`
    fastq_reads_tmp=`cat {fastq_line_num_file}`
    fastq_reads=`expr $fastq_reads_tmp / 2`

    if [ $total_reads -ne $fastq_reads ]; then
        echo "Total read count is not good for this data. BAM file: ${{total_reads}} reads. FASTQ file: ${{fastq_reads}} reads." >&2
        exit 1
    fi
fi

genomon_qc merge {coverage_file} {bamstats_file} {output_file} --meta "{meta}"
"""

def configure(input_bams, genomon_conf, run_conf, sample_conf):
    
    STAGE_NAME = "qc_merge"
    CONF_SECTION = "qc_merge"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": genomon_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": genomon_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Qc_merge(params)
    
    output_files = []
    for sample in sample_conf.qc:
        output_file = "qc/%s/%s.genomonQC.result.txt" % (sample, sample)
        output_files.append(output_file)
        
        arguments = {
            "bamstats_file": "%s/qc/%s/%s.bamstats" % (run_conf.project_root, sample, sample),
            "coverage_file": "%s/qc/%s/%s.coverage" % (run_conf.project_root, sample, sample),
            "output_file": "%s/%s" % (run_conf.project_root, output_file),
            "meta": genomon_conf.get_meta_info(["genomon_pipeline", "qc_coverage"]),
            "fastq_line_num_file": run_conf.project_root +'/fastq/'+ sample +'/fastq_line_num.txt'
        }

        singularity_bind = [run_conf.project_root]

        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
