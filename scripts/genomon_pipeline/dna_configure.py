import os

## link the input fastq to project directory
def link_input_fastq(sample, output_file, genomon_conf, run_conf, sample_conf):
    fastq_dir = run_conf.project_root + '/fastq/' + sample
    fastq_prefix, ext = os.path.splitext(sample_conf.fastq[sample][0][0])
    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    for (count, fastq_files) in enumerate(sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_files)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_1'+ ext):
            os.symlink(sample_conf.fastq[sample][0][count], fastq_dir + '/'+str(count+1)+'_1'+ ext)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_2'+ ext):
            os.symlink(sample_conf.fastq[sample][1][count], fastq_dir + '/'+str(count+1)+'_2'+ ext)

# merge sorted bams into one and mark duplicate reads with biobambam
def bwa(input_files, output_file, output_dir, genomon_conf, run_conf, sample_conf):
    import genomon_pipeline.dna_resource.markduplicates  as dr_markduplicates 
    markduplicates = dr_markduplicates.Markduplicates(genomon_conf.get("markduplicates", "qsub_option"), run_conf.drmaa)

    sample_name = os.path.basename(output_dir)

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file}
    
    singularity_params = {
        "image": genomon_conf.get("markduplicates", "image"),
        "option": genomon_conf.get("markduplicates", "singularity_option"),
        "bind": [run_conf.project_root],
    }
        
    markduplicates.task_exec(
            arguments, run_conf.project_root + '/log/' + sample_name , 
            run_conf.project_root + '/script/'+ sample_name, singularity_params)

    for input_file in input_files:
        os.unlink(input_file)
        os.unlink(input_file + ".bai")

def main(genomon_conf, run_conf, sample_conf):
    
    # generate output list of 'linked fastq'
    linked_fastq_list = []
    for sample in sample_conf.fastq:
        if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
        if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue
    
        link_fastq_arr1 = []
        link_fastq_arr2 = []
        for (count, fastq_file) in enumerate(sample_conf.fastq[sample][0]):
            fastq_prefix, ext = os.path.splitext(fastq_file)
            link_fastq_arr1.append(run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_1' + ext)
            link_fastq_arr2.append(run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_2' + ext)
        #linked_fastq_list.append([link_fastq_arr1,link_fastq_arr2])
        link_input_fastq(sample, [link_fastq_arr1,link_fastq_arr2], genomon_conf, run_conf, sample_conf)
        
    # generate output list of 'bam2fastq'
    bam2fastq_output_list = []
    for sample in sample_conf.bam_tofastq:
        if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
        if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue
        bam2fastq_arr1 = []
        bam2fastq_arr2 = []
        bam2fastq_arr1.append(run_conf.project_root + '/fastq/' + sample + '/1_1.fastq')
        bam2fastq_arr2.append(run_conf.project_root + '/fastq/' + sample + '/1_2.fastq')
        bam2fastq_output_list.append([bam2fastq_arr1,bam2fastq_arr2])
    
   