#! /usr/bin/env python

import os
import yaml

class Stage_task(object):
    params = {
        "work_dir": "",
        "stage_name": "",
        "image": "",
        "qsub_option": "",
        "singularity_option": ""
    }
    
    def __init__(self, params):
        self.script_dir = "%s/script" % (os.path.abspath(params["work_dir"]))
        self.shell_script_name = "shell_%s.sh" % (params["stage_name"])
        self.singurality_script_name = "singularity_%s.sh" % (params["stage_name"])
        self.snakemake_conf_name = "conf_%s.yml" % (params["stage_name"])
        
        self.image = params["image"]
        self.qsub_option = params["qsub_option"]
        self.singularity_option = params["singularity_option"]
        self.log_dir = "%s/log" % (params["work_dir"])
        
        self.shell_script_template = ""
        self.singurality_script_template = """#!/bin/bash
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

singularity exec {option} --bind {bind} {image} /bin/bash {script}
"""
        
    def write_script(self, arguments, singularity_bind, run_conf, sample = "", max_task = 0):
        output_dir = self.script_dir
        log_dir = self.log_dir
        if sample != "":
            output_dir = "%s/%s" % (self.script_dir, sample)
            log_dir  = "%s/%s" % (self.log_dir, sample)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(log_dir, exist_ok=True)

        shell_script_path = "%s/%s" % (output_dir, self.shell_script_name)
        open(shell_script_path, 'w').write(self.shell_script_template.format(**arguments))
        
        singurality_script_path = "%s/%s" % (output_dir, self.singurality_script_name)
        open(singurality_script_path, "w").write(self.singurality_script_template.format(
            option = self.singularity_option,
            bind = ",".join(singularity_bind),
            image = self.image,
            script = shell_script_path
        ))
        
        conf_path = "%s/%s" % (output_dir, self.snakemake_conf_name)
        open(conf_path, "w").write(yaml.dump({
            "qsub_option": self.qsub_option,
            "log_dir": log_dir,
            "drmaa": run_conf.drmaa,
            "retry_count": run_conf.retry_count,
            "max_task": max_task
        }))
        
    def configure(self):
        pass
