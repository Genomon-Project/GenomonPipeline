#! /usr/bin/env python

import genomon_pipeline.core.stage_task_abc as stage_task

class Qc_coverage(stage_task.Stage_task):
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
set -o pipefail

# set python environment
export bedtools=/usr/local/bin/bedtools
export samtools=/usr/local/bin/samtools

if [ "{data_type}" = "wgs" ]
then
########## WGS ##########
genomon_qc wgs {input_file} {output_file} --genome_size_file {genome_size_file} --gaptxt {gaptxt} --incl_bed_width {incl_bed_width} --i_bed_lines {i_bed_lines} --i_bed_width {i_bed_width} --bedtools $bedtools --samtools $samtools --samtools_params "{samtools_params}" --coverage_text {coverage_text} {grc_flag} --del_tempfile

else
########## exome ##########
genomon_qc exome {input_file} {output_file} --bait_file {bait_file} --bedtools $bedtools --samtools $samtools --samtools_params "{samtools_params}" --coverage_text {coverage_text} {grc_flag} --del_tempfile
fi
"""

def configure(input_bams, genomon_conf, run_conf, sample_conf):
    
    STAGE_NAME = "qc_coverage"
    CONF_SECTION = "qc_coverage"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": genomon_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": genomon_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Qc_coverage(params)
    
    data_type = "exome"
    if genomon_conf.get("qc_coverage", "wgs_flag") == "True":
        data_type = "wgs"
    
    grc_flag = ""
    if genomon_conf.get("qc_coverage", "grc_flag") == "True":
        grc_flag = "--grc_flag"

    output_files = []
    for sample in sample_conf.qc:
        output_file = "qc/%s/%s.coverage" % (sample, sample)
        output_files.append(output_file)

        arguments = {
            "data_type": data_type,
            "coverage_text": genomon_conf.get("qc_coverage", "coverage"),
            "i_bed_lines": genomon_conf.get("qc_coverage", "wgs_i_bed_lines"),
            "i_bed_width": genomon_conf.get("qc_coverage", "wgs_i_bed_width"),
            "incl_bed_width": genomon_conf.get("qc_coverage", "wgs_incl_bed_width"),
            "genome_size_file": genomon_conf.path_get("qc_coverage", "genome_size"),
            "gaptxt": genomon_conf.path_get("qc_coverage", "gaptxt"),
            "bait_file": genomon_conf.path_get("qc_coverage", "bait_file"),
            "samtools_params": genomon_conf.get("qc_coverage", "samtools_params"),
            "grc_flag": grc_flag,
            "input_file": input_bams[sample],
            "output_file": "%s/%s" % (run_conf.project_root, output_file)
        }

        singularity_bind = [
            run_conf.project_root,
            genomon_conf.path_get("qc_coverage", "genome_size"),
            genomon_conf.path_get("qc_coverage", "gaptxt"),
            genomon_conf.path_get("qc_coverage", "bait_file")
        ]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]

        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
