#! /usr/bin/env python

def main(genomon_conf, run_conf, sample_conf):
    
    # preparation
    import genomon_pipeline.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq_single, sample_conf.bam_tofastq_pair)
    setup.create_directories(genomon_conf, run_conf, input_stages, 'rna/data/snakefile.txt')
    bam_tofastq_stages = (sample_conf.bam_tofastq_single, sample_conf.bam_tofastq_pair)
    setup.touch_bam_tofastq(run_conf, bam_tofastq_stages)
    
    # dump conf.yaml
    import genomon_pipeline.rna.resource.star_align as rs_align
    y = setup.dump_yaml_input_section(
        run_conf,
        (sample_conf.bam_tofastq_single, sample_conf.bam_tofastq_pair),
        sample_conf.fastq,
        sample_conf.bam_import, 
        rs_align.OUTPUT_FORMAT
    )

    # link fastq
    linked_fastqs = setup.link_input_fastq(run_conf, sample_conf.fastq, sample_conf.fastq_src)
    
    # link import bam
    output_bams = setup.link_import_bam(
        run_conf, sample_conf.bam_import, 
        rs_align.BAM_POSTFIX, 
        rs_align.BAI_POSTFIX,
        rs_align.OUTPUT_FORMAT.split("/")[0]
    )
    
    # ######################
    # create scripts
    # ######################
    # bam to fastq
    import genomon_pipeline.rna.resource.bamtofastq_single as rs_bamtofastq_single
    output_fastqs = rs_bamtofastq_single.configure(genomon_conf, run_conf, sample_conf)
    import genomon_pipeline.rna.resource.bamtofastq_pair as rs_bamtofastq_pair
    output_fastqs.update(rs_bamtofastq_pair.configure(genomon_conf, run_conf, sample_conf))
    
    # star
    for sample in output_fastqs:
        sample_conf.fastq[sample] = output_fastqs[sample]
        sample_conf.fastq_src[sample] = []

    for sample in linked_fastqs:
        sample_conf.fastq[sample] = linked_fastqs[sample]["fastq"]
        sample_conf.fastq_src[sample] = linked_fastqs[sample]["src"]

    align_bams = rs_align.configure(genomon_conf, run_conf, sample_conf)
    output_bams.update(align_bams)

    # expression
    import genomon_pipeline.rna.resource.expression as rs_expression
    output_expressions = rs_expression.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    # ######################
    # dump conf.yaml
    # ######################
    y["output_files"].extend(output_expressions)
    
    y["expression_samples"] = {}
    for sample in sample_conf.expression:
        y["expression_samples"][sample] = rs_align.OUTPUT_FORMAT.format(sample=sample)
        
    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))
    