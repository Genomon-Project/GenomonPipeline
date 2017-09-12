#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_QC_Coverage(Stage_task):

    task_name = "qc_coverage"

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
set -eu
set -o pipefail

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export PYTHONPATH={pythonpath}

if [ {data_type} = "wgs" ]
then
########## WGS ##########
{genomon_qc} wgs {input_file} {output_file} --genome_size_file {genome_size_file} --gaptxt {gaptxt} --incl_bed_width {incl_bed_width} --i_bed_lines {i_bed_lines} --i_bed_width {i_bed_width} --ld_library_path {ld_library_path} --bedtools {bedtools} --samtools {samtools} --samtools_params "{samtools_params}" --coverage_text {coverage_text} --del_tempfile True

else
########## exome ##########
{genomon_qc} exome {input_file} {output_file} --bait_file {bait_file} --ld_library_path {ld_library_path} --bedtools {bedtools} --samtools {samtools} --samtools_params "{samtools_params}" --coverage_text {coverage_text} --del_tempfile True
fi
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Coverage, self).__init__(qsub_option, script_dir)
