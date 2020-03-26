#! /usr/bin/env python

import genomon_pipeline.core.stage_task_abc as stage_task

class Bam_tofastq_single(stage_task.Stage_task):
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
export LD_LIBRARY_PATH=/usr/local/lib
bams=( `echo "{input_bam}" | tr -s ';' ' '`)
if [ ${{#bams[@]}} -eq 1 ]; then
    bam=${{bams[0]}}
    echo $bam
    /usr/local/bin/bamtofastq {param} filename=${{bam}} F={f1_name} T={t} S={s} O={o1_name} || exit $?    
else
    echo -n > {f1_name}
    echo -n > {t}
    echo -n > {s}
    echo -n > {o1_name}
    for bam in ${{bams[@]}}; do
        echo $bam
        /usr/local/bin/bamtofastq {param} filename=${{bam}} F={f1_name}.tmp T={t}.tmp S={s}.tmp O={o1_name}.tmp || exit $?
        cat {f1_name}.tmp >> {f1_name} || exit $?
        if [ -s {t}.tmp ]; then
            cat {t}.tmp >> {t} || exit $?
        fi
        if [ -s {s}.tmp ]; then
            cat {s}.tmp >> {s} || exit $?
        fi
        if [ -s {o1_name}.tmp ]; then
            cat {o1_name}.tmp >> {o1_name} || exit $?
        fi
        rm {f1_name}.tmp
        rm {t}.tmp
        rm {s}.tmp
        rm {o1_name}.tmp
    done
fi
touch {pass}
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(genomon_conf, run_conf, sample_conf):
    
    #STAGE_NAME = sample_conf.SECTION_BAM_TOFASTQ_SINGLE
    STAGE_NAME = "bam_tofastq"
    CONF_SECTION = "bam_tofastq"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": genomon_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": genomon_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": genomon_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Bam_tofastq_single(params)
    
    output_fastqs = {}
    for sample in sample_conf.bam_tofastq_single:
        output_dir = "%s/fastq/%s" % (run_conf.project_root, sample)
        f1_name = output_dir + "/1_1.fastq"
        output_fastqs[sample] = [[f1_name]]

        arguments = {
            "param": genomon_conf.get(CONF_SECTION, "params"),
            "input_bam": sample_conf.bam_tofastq_single[sample],
            "f1_name": f1_name,
            "o1_name": output_dir + '/unmatched_first_output.txt',
            "t": output_dir + '/temp.txt',
            "s": output_dir + '/single_end_output.txt',
            "pass": output_dir + '/pass.txt'
        }
        
        singularity_bind = [
            run_conf.project_root,
        ] + sample_conf.bam_tofastq_single_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_fastqs
