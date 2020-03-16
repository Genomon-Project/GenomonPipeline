#! /usr/bin/env python

import os
import genomon_pipeline.config.genomon_conf as gc
import genomon_pipeline.config.run_conf as rc
import genomon_pipeline.config.sample_conf as sc

def main(args):

    ###
    # set run_conf
    run_conf = rc.Run_conf()
    run_conf.sample_conf_file = args.sample_conf_file
    run_conf.analysis_type = args.analysis_type
    run_conf.project_root = os.path.abspath(args.project_root)
    run_conf.genomon_conf_file = args.genomon_conf_file
    run_conf.drmaa = False if args.disable_drmaa else True
    run_conf.retry_count = args.retry_count
    
    ###
    # read sample list file
    sample_conf = sc.Sample_conf(run_conf.sample_conf_file)
    
    ###
    # set genomon_conf and task parameter config data
    genomon_conf = gc.Genomon_conf(conf = run_conf.genomon_conf_file)
    
    if run_conf.analysis_type == "dna":
        #genomon_conf.genomon_conf_check()
        genomon_conf.software_version_set()
        import genomon_pipeline.dna_configure
        genomon_pipeline.dna_configure.main(genomon_conf = genomon_conf, run_conf = run_conf, sample_conf = sample_conf)
        