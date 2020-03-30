#! /usr/bin/env python

import genomon_pipeline.core.stage_task_abc as stage_task

class Qc_bamstats(stage_task.Stage_task):
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
set -o pipefail

# set python environment
export bamstats=/tools/ICGC/bin/bam_stats.pl
export ld_library_path=/usr/local/bin
export perl5lib=/tools/ICGC/lib/perl5:/tools/ICGC/lib/perl5/x86_64-linux-gnu-thread-multi

mkdir -p $(dirname {OUTPUT_FILE})
genomon_qc bamstats {INPUT_BAM} {OUTPUT_FILE} --perl5lib $perl5lib --bamstats $bamstats
"""

def configure(input_bams, genomon_conf, run_conf, sample_conf):
    
    STAGE_NAME = "qc_bamstats"
    CONF_SECTION = "qc_bamstats"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": genomon_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": genomon_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Qc_bamstats(params)
    
    output_files = []
    for sample in sample_conf.qc:
        output_file = "qc/%s/%s.bamstats" % (sample, sample)
        output_files.append(output_file)
        
        arguments = {
            "INPUT_BAM": input_bams[sample],
            "OUTPUT_FILE": "%s/%s" % (run_conf.project_root, output_file),
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]

        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
