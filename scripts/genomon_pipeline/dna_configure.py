import os
import shutil
import pkg_resources

def create_directories(genomon_conf, run_conf, sample_conf):
    os.makedirs(run_conf.project_root + '/config/', exist_ok=True)
    genomon_conf_name, genomon_conf_ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
    sample_conf_name, sample_conf_ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))
    shutil.copyfile(run_conf.genomon_conf_file, run_conf.project_root + '/config/' + genomon_conf_name +'_'+ genomon_conf.analysis_timestamp + genomon_conf_ext)
    shutil.copyfile(run_conf.sample_conf_file, run_conf.project_root + '/config/' + sample_conf_name +'_'+ genomon_conf.analysis_timestamp + sample_conf_ext)
    shutil.copyfile(pkg_resources.resource_filename('genomon_pipeline', 'data/snakefile'), run_conf.project_root + '/snakefile')
    
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
            ext = ".fastq"
            if not os.path.exists(fastq_dir + '/'+str(count+1)+'_1'+ ext):
                os.symlink(sample_conf.fastq[sample][0][count], fastq_dir + '/'+str(count+1)+'_1'+ ext)
            if not os.path.exists(fastq_dir + '/'+str(count+1)+'_2'+ ext):
                os.symlink(sample_conf.fastq[sample][1][count], fastq_dir + '/'+str(count+1)+'_2'+ ext)

# link the import bam to project directory
def link_import_bam(genomon_conf, run_conf, sample_conf):
    linked_bam = {}
    for sample in sample_conf.bam_import:
        bam = sample_conf.bam_import[sample]
        link_dir = run_conf.project_root + '/bam/' + sample
        os.makedirs(link_dir, exist_ok=True)
        bam_prefix, ext = os.path.splitext(bam)
        linked_bam[sample] = link_dir +'/'+ sample +'.markdup.bam'

        if (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam')) and (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam.bai')): 
            os.symlink(bam, link_dir +'/'+ sample +'.markdup.bam')
            if (os.path.exists(bam +'.bai')):
                os.symlink(bam +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')
            elif (os.path.exists(bam_prefix +'.bai')):
                os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')
    return linked_bam

def dump_conf_yaml(genomon_conf, run_conf, sample_conf):
    samples = []
    outputs =[]
    for sample in sample_conf.fastq:
        samples.append(sample)
        outputs.append("bam/{sample}/{sample}.markdup.bam".format(sample = sample))
    
    input_mutation = {}
    for (sample, control, control_panel) in sample_conf.mutation_call:
        input_mutation[sample] = ("bam/%s/%s.markdup.bam" % (sample, sample))
        outputs.append("mutation/%s/%s.txt" % (sample, sample))
    
    input_sv = {}
    for (sample, control, control_panel) in sample_conf.sv_detection:
        input_sv[sample] = ("bam/%s/%s.markdup.bam" % (sample, sample))
        outputs.append("sv/%s/%s.txt" % (sample, sample))
    
    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump({
        "samples": samples,
        "mutation_samples": input_mutation,
        "sv_samples": input_sv,
        "output_files": outputs
    }))
    
def main(genomon_conf, run_conf, sample_conf):
 
    create_directories(genomon_conf, run_conf, sample_conf)
    link_input_fastq(genomon_conf, run_conf, sample_conf)
    output_bams = link_import_bam(genomon_conf, run_conf, sample_conf)
    
    import genomon_pipeline.resource.bwa_align
    align_bams = genomon_pipeline.resource.bwa_align.configure(genomon_conf, run_conf, sample_conf)
    output_bams.update(align_bams)
    
    import genomon_pipeline.resource.mutation_dummy
    genomon_pipeline.resource.mutation_dummy.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    import genomon_pipeline.resource.sv_dummy
    genomon_pipeline.resource.sv_dummy.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    dump_conf_yaml(genomon_conf, run_conf, sample_conf)
