#! /usr/bin/env python

import os
import shutil
import pkg_resources

def create_directories(genomon_conf, run_conf, input_stages, snakefile_name):
    if not type(input_stages) in [type(()), type([])]:
        raise Exception("type of input_stages must be tuple or list")
        
    # mkdir config
    os.makedirs(run_conf.project_root + '/config/', exist_ok=True)
    
    # copy genomon.cfg
    genomon_conf_name, genomon_conf_ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
    shutil.copyfile(run_conf.genomon_conf_file, run_conf.project_root + '/config/' + genomon_conf_name +'_'+ genomon_conf.analysis_timestamp + genomon_conf_ext)
    
    # copy sample.csv
    sample_conf_name, sample_conf_ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))
    shutil.copyfile(run_conf.sample_conf_file, run_conf.project_root + '/config/' + sample_conf_name +'_'+ genomon_conf.analysis_timestamp + sample_conf_ext)
    
    # copy snakemake
    shutil.copyfile(pkg_resources.resource_filename('genomon_pipeline', snakefile_name), run_conf.project_root + '/snakefile')
    
    # mkdir log
    for stage in input_stages:
        for sample in stage:
            os.makedirs(run_conf.project_root + '/log/' + sample, exist_ok=True)

# touch snakemake entry-file
def touch_bam_tofastq(genomon_conf, run_conf, bam_tofastq_stages):
    if not type(bam_tofastq_stages) in [type(()), type([])]:
        raise Exception("type of bam_tofastq_stages must be tuple or list")
        
    for stage in bam_tofastq_stages:
        for sample in stage:
            wdir = run_conf.project_root + '/bam_tofastq/' + sample
            os.makedirs(wdir, exist_ok=True)
            open(wdir + '/' + sample + ".txt", "w").close()

# link the input fastq to project directory
def link_input_fastq(genomon_conf, run_conf, fastq_stage, fastq_stage_src):
    linked_fastq = {}
    for sample in fastq_stage:
        fastq_dir = run_conf.project_root + '/fastq/' + sample
        os.makedirs(fastq_dir, exist_ok=True)
        open("%s/pass.txt" % (fastq_dir), "w").close()
        pair = len(fastq_stage[sample]) > 1
            
        new_fastq_src = []
        new_fastq_src += fastq_stage_src[sample]
        new_fastq_src += fastq_stage[sample][0]
        if pair:
            new_fastq_src += fastq_stage[sample][1]

        linked_fastq[sample] = {"fastq": [[], []], "src": new_fastq_src}
        
        for (count, fastq_files) in enumerate(fastq_stage[sample][0]):
            fastq_prefix, ext = os.path.splitext(fastq_files)
            r1 = fastq_dir + '/'+str(count+1)+'_1'+ ext
            linked_fastq[sample]["fastq"][0] += [r1]
            if not os.path.exists(r1):
                os.symlink(fastq_stage[sample][0][count], r1)
            
            if pair:
                r2 = fastq_dir + '/'+str(count+1)+'_2'+ ext
                linked_fastq[sample]["fastq"][1] += [r2]
                if not os.path.exists(r2):
                    os.symlink(fastq_stage[sample][1][count], r2)

    return linked_fastq

# link the import bam to project directory
def link_import_bam(genomon_conf, run_conf, bam_import_stage, bam_postfix, bai_postfix, subdir = "bam"):
    linked_bam = {}
    for sample in bam_import_stage:
        bam = bam_import_stage[sample]
        link_dir = "%s/%s/%s" % (run_conf.project_root, subdir, sample)
        os.makedirs(link_dir, exist_ok=True)
        prefix, ext = os.path.splitext(bam)
        if ext == ".bam":
            bai = ".bai"
        elif ext == ".cram":
            bai = ".crai"
        linked_bam[sample] = link_dir +'/'+ sample + bam_postfix

        if (not os.path.exists(link_dir +'/'+ sample + bam_postfix)) and (not os.path.exists(link_dir +'/'+ sample + bai_postfix)): 
            os.symlink(bam, link_dir +'/'+ sample + bam_postfix)
            if (os.path.exists(bam + bai)):
                os.symlink(bam + bai, link_dir +'/'+ sample + bai_postfix)
            elif (os.path.exists(bam_postfix + bai)):
                os.symlink(bam_postfix + bai, link_dir +'/'+ sample + bai_postfix)
    return linked_bam


def dump_yaml_input_section(genomon_conf, run_conf, bam_tofastq_stages, fastq_stage, bam_import_stage, bam_template):
    if not type(bam_tofastq_stages) in [type(()), type([])]:
        raise Exception("type of bam_tofastq_stages must be tuple or list")
    
    samples = []
    outputs =[]
    for stage in bam_tofastq_stages:
        for sample in stage:
            samples.append(sample)
            outputs.append("fastq/%s/pass.txt" % (sample))
        
    input_aln = {}
    for sample in fastq_stage:
        input_aln[sample] = "fastq/%s/pass.txt" % (sample)
        outputs.append(bam_template.format(sample = sample))
    
    for sample in bam_import_stage:
        input_aln[sample] = "fastq/%s/pass.txt" % (sample)
        
    return {
        "samples": samples,
        "aln_samples": input_aln,
        "output_files": outputs
    }
