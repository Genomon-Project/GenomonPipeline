#! /usr/bin/env python

import os
import genomon_pipeline.config.genomon_conf as gc
import genomon_pipeline.config.run_conf as rc
import genomon_pipeline.config.sample_conf as sc

def main(args):

    ###
    # set run_conf
    rc.run_conf.sample_conf_file = args.sample_conf_file
    rc.run_conf.analysis_type = args.analysis_type
    rc.run_conf.project_root = os.path.abspath(args.project_root)
    rc.run_conf.genomon_conf_file = args.genomon_conf_file
    rc.run_conf.drmaa = False if args.disable_drmaa else True

    ###
    # read sample list file
    sc.sample_conf.parse_file(rc.run_conf.sample_conf_file)

    ###
    # set genomon_conf and task parameter config data
    gc.genomon_conf.read(rc.run_conf.genomon_conf_file)
    
    if rc.run_conf.analysis_type == "dna":
        gc.dna_genomon_conf_check()
        gc.dna_software_version_set()
        import genomon_pipeline.dna_pipeline
        genomon_pipeline.dna_pipeline.configure(genomon_conf = gc, run_conf = rc, sample_conf = sc)
        