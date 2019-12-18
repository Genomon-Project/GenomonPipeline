#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Res_PostAnalysis(st.Stage_task):

    task_name = "post_analysis"

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

samtools=samtools
bedtools=bedtools

genomon_pa {mode} {output_dir} {genomon_root} {sample_sheet} \
--config_file {config_file} \
--samtools $samtools --bedtools $bedtools \
--input_file_case1 "{input_file_case1}" \
--input_file_case2 "{input_file_case2}" \
--input_file_case3 "{input_file_case3}" \
--input_file_case4 "{input_file_case4}"
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PostAnalysis, self).__init__(qsub_option, script_dir)

    def output_files(self, mode, samples, genomon_root, sample_conf_name, genomon_conf):
        
        import os
        import sys
        if sys.version_info.major == 2:
            import ConfigParser as cp
        else:
            import configparser as cp
        
        pa_conf = cp.RawConfigParser()
        pa_conf.read(genomon_conf.get("post_analysis", "config_file"))

        analysis_dir = genomon_root + "/" + mode
        pa_dir = genomon_root + "/post_analysis/" + sample_conf_name
        
        if mode == "qc":
            di_outputs = { "outputs": [], "run_pa":False}
        
            if len(samples) == 0:
                return di_outputs
            
            run_analysis = False
            
            for sample in samples:
                if not os.path.exists(analysis_dir + "/" + sample + "/" + sample + pa_conf.get("result_format_qc", "suffix")):
                    run_analysis = True
                    break
            
            if run_analysis == True: di_outputs["run_pa"] = True
            
            output = pa_dir + '/' + pa_conf.get("merge_format_qc", "filename_all")
            if not os.path.exists(output): di_outputs["run_pa"] = True
            
            di_outputs["outputs"].append(output)
            
            return di_outputs
            
        section = ""
        section_in = ""
        
        if mode == "mutation":
            section = "merge_format_mutation"
            section_in = "result_format_mutation"
        elif mode == "sv":
            section = "merge_format_sv"
            section_in = "result_format_sv"
        else:
            return {}
            
        di_outputs = { \
            "case1":{"output_filt":"", "output_unfilt":"", "run_pa":False, "samples":[]}, \
            "case2":{"output_filt":"", "output_unfilt":"", "run_pa":False, "samples":[]}, \
            "case3":{"output_filt":"", "output_unfilt":"", "run_pa":False, "samples":[]}, \
            "case4":{"output_filt":"", "output_unfilt":"", "run_pa":False, "samples":[]}, \
            "all":  {"output_filt":"", "output_unfilt":"", "run_pa":False}, \
            "outputs": [], \
            "run_pa": False \
        }
        
        if len(samples) == 0:
            return di_outputs

        run_analysis = {"case1":False, "case2":False, "case3":False, "case4":False, "all":  False}
        
        use_case1 = pa_conf.getboolean(section, "output_case1") or pa_conf.getboolean(section, "output_filt_case1")
        use_case2 = pa_conf.getboolean(section, "output_case2") or pa_conf.getboolean(section, "output_filt_case2")
        use_case3 = pa_conf.getboolean(section, "output_case3") or pa_conf.getboolean(section, "output_filt_case3")
        use_case4 = pa_conf.getboolean(section, "output_case4") or pa_conf.getboolean(section, "output_filt_case4")
        output_all = pa_conf.getboolean(section, "output_all") or pa_conf.getboolean(section, "output_filt_all")
        
        for complist in samples:
            type = ""
            if   (complist[1] != None and complist[2] != None and use_case1 == True): type = "case1"
            elif (complist[1] != None and complist[2] == None and use_case2 == True): type = "case2"
            elif (complist[1] == None and complist[2] != None and use_case3 == True): type = "case3"
            elif (complist[1] == None and complist[2] == None and use_case4 == True): type = "case4"
            else: continue
            
            di_outputs[type]["samples"].append(complist[0])
            if not os.path.exists(analysis_dir + "/" + complist[0] + "/" + complist[0] + pa_conf.get(section_in, "suffix")):
                # print analysis_dir + "/" + complist[0] + "/" + complist[0] + pa_conf.get(section_in, "suffix")
                run_analysis[type] = True
                if output_all:
                    run_analysis["all"] = True
        
        # each case
        for key in ["case1", "case2", "case3", "case4", "all"]:
            if key != "all":
                if len(di_outputs[key]["samples"]) == 0: continue
            
            if run_analysis[key] == True:
                di_outputs[key]["run_pa"] = True
            
            if pa_conf.getboolean(section, "output_filt_" + key):
                output_filt = pa_dir + '/' + pa_conf.get(section, "filename_filt_" + key)
                di_outputs[key]["output_filt"] = output_filt
                di_outputs["outputs"].append(output_filt)
                
                if not os.path.exists(output_filt):
                    di_outputs[key]["run_pa"] = True
            
            if pa_conf.getboolean(section, "output_" + key):
                output_unfilt = pa_dir + '/' + pa_conf.get(section, "filename_" + key)
                di_outputs[key]["output_unfilt"] = output_unfilt
                di_outputs["outputs"].append(output_unfilt)
                
                if not os.path.exists(output_unfilt):
                    di_outputs[key]["run_pa"] = True
            
            if di_outputs[key]["run_pa"] == True:
                di_outputs["run_pa"] = True
            
        return di_outputs
