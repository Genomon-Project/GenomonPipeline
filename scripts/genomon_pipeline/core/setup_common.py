#! /usr/bin/env python

import os
import shutil
import pkg_resources

def create_directories(genomon_conf, run_conf, sample_conf, snakefile_name):
    # mkdir config
    os.makedirs(run_conf.project_root + '/config/', exist_ok=True)
    
    # copy genomon.cfg
    genomon_conf_name, genomon_conf_ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
    shutil.copyfile(run_conf.genomon_conf_file, run_conf.project_root + '/config/' + genomon_conf_name +'_'+ genomon_conf.analysis_timestamp + genomon_conf_ext)
    
    # copy sample.csv
    sample_conf_name, sample_conf_ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))
    shutil.copyfile(run_conf.sample_conf_file, run_conf.project_root + '/config/' + sample_conf_name +'_'+ genomon_conf.analysis_timestamp + sample_conf_ext)
    
    # copy snakemake
    shutil.copyfile(pkg_resources.resource_filename('genomon_pipeline', 'data/' + snakefile_name), run_conf.project_root + '/snakefile')
    
    # mkdir loag
    for target_sample_dict in (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq):
        for sample in target_sample_dict:
            os.makedirs(run_conf.project_root + '/log/' + sample, exist_ok=True)
    
# link the input fastq to project directory
def link_input_fastq(genomon_conf, run_conf, sample_conf):
    for sample in sample_conf.fastq:
        fastq_dir = run_conf.project_root + '/fastq/' + sample
        os.makedirs(run_conf.project_root + '/fastq/' + sample, exist_ok=True)
        
        for (count, fastq_files) in enumerate(sample_conf.fastq[sample][0]):
            fastq_prefix, ext = os.path.splitext(fastq_files)
            #ext = ".fastq"
            if not os.path.exists(fastq_dir + '/'+str(count+1)+'_1'+ ext):
                os.symlink(sample_conf.fastq[sample][0][count], fastq_dir + '/'+str(count+1)+'_1'+ ext)
            if not os.path.exists(fastq_dir + '/'+str(count+1)+'_2'+ ext):
                os.symlink(sample_conf.fastq[sample][1][count], fastq_dir + '/'+str(count+1)+'_2'+ ext)
            
            open(fastq_dir + '/' + sample + ".txt", "w").close()
                
# link the import bam to project directory
def link_import_bam(genomon_conf, run_conf, sample_conf, bam_prefix, bai_prefix):
    linked_bam = {}
    for sample in sample_conf.bam_import:
        bam = sample_conf.bam_import[sample]
        link_dir = run_conf.project_root + '/bam/' + sample
        os.makedirs(link_dir, exist_ok=True)
        bam_prefix, ext = os.path.splitext(bam)
        linked_bam[sample] = link_dir +'/'+ sample + bam_prefix

        if (not os.path.exists(link_dir +'/'+ sample + bam_prefix)) and (not os.path.exists(link_dir +'/'+ sample + bai_prefix)): 
            os.symlink(bam, link_dir +'/'+ sample + bam_prefix)
            if (os.path.exists(bam +'.bai')):
                os.symlink(bam +'.bai', link_dir +'/'+ sample + bai_prefix)
            elif (os.path.exists(bam_prefix +'.bai')):
                os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample + bai_prefix)
    return linked_bam

