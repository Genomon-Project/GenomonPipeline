#! /usr/bin/env python

import os
import genomon_pipeline.core.genomon_conf as gc
import genomon_pipeline.core.run_conf as rc

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

    # disable params
    # - args.multiprocess
    # - args.param_check
    
    ###
    # set genomon_conf and task parameter config data
    genomon_conf = gc.Genomon_conf(conf = run_conf.genomon_conf_file)
    genomon_conf.software_version_set()
    
    if run_conf.analysis_type == "dna":
        import genomon_pipeline.dna.sample_conf as sc
        import genomon_pipeline.dna.configure as configure

    elif run_conf.analysis_type == "rna":
        import genomon_pipeline.rna.sample_conf as sc
        import genomon_pipeline.rna.configure as configure
    
    elif run_conf.analysis_type == "germ":
        import genomon_pipeline.germ.sample_conf as sc
        import genomon_pipeline.germ.configure as configure
        
    sample_conf = sc.Sample_conf(run_conf.sample_conf_file)
    configure.main(genomon_conf = genomon_conf, run_conf = run_conf, sample_conf = sample_conf)

