#! /usr/bin/env python

def main(genomon_conf, run_conf, sample_conf):
    
    # preparation
    import genomon_pipeline.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq)
    setup.create_directories(genomon_conf, run_conf, input_stages, 'germ/data/snakefile.txt')
    setup.touch_bam_tofastq(run_conf, (sample_conf.bam_tofastq,))
    
    # dump conf.yaml
    import genomon_pipeline.germ.resource.fq2cram as rs_align
    y = setup.dump_yaml_input_section(
        run_conf, 
        (sample_conf.bam_tofastq, ),
        sample_conf.fastq,
        sample_conf.bam_import, 
        rs_align.OUTPUT_FORMAT
    )

    # link fastq
    linked_fastqs = setup.link_input_fastq(run_conf, sample_conf.fastq, sample_conf.fastq_src)
    
    # link import bam
    output_bams = setup.link_import_bam(
        run_conf,
        sample_conf.bam_import, 
        rs_align.BAM_POSTFIX, 
        rs_align.BAI_POSTFIX,
        rs_align.OUTPUT_FORMAT.split("/")[0]
    )
    
    # ######################
    # create scripts
    # ######################
    # bam to fastq
    import genomon_pipeline.germ.resource.bamtofastq as rs_bamtofastq
    output_fastqs = rs_bamtofastq.configure(genomon_conf, run_conf, sample_conf)
    
    # bwa
    for sample in output_fastqs:
        sample_conf.fastq[sample] = output_fastqs[sample]
        sample_conf.fastq_src[sample] = []

    for sample in linked_fastqs:
        sample_conf.fastq[sample] = linked_fastqs[sample]["fastq"]
        sample_conf.fastq_src[sample] = linked_fastqs[sample]["src"]
    
    align_bams = rs_align.configure(genomon_conf, run_conf, sample_conf)
    output_bams.update(align_bams)

    # mutation
    import genomon_pipeline.germ.resource.haplotypecaller as rs_mutation
    output_mutations = rs_mutation.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    
    # ######################
    # dump conf.yaml
    # ######################
    y["output_files"].extend(output_mutations)
    
    y["htc_samples"] = {}
    for sample in sample_conf.haplotype_call:
        y["htc_samples"][sample] = rs_align.OUTPUT_FORMAT.format(sample=sample)
        
    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))
    