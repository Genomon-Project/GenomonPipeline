#! /usr/bin/env python

def dump_conf_yaml(genomon_conf, run_conf, sample_conf):
    
    import genomon_pipeline.core.setup_common as setup
    
    y = setup.dump_yaml_input_section(
        genomon_conf, 
        run_conf,
        (sample_conf.bam_tofastq_single, sample_conf.bam_tofastq_pair),
        sample_conf.fastq,
        sample_conf.bam_import, 
        "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    )
    input_expression = {}
    for sample in sample_conf.expression:
        input_expression[sample] = "star/{sample}/{sample}.Aligned.sortedByCoord.out.bam".format(sample = sample)
        y["output_files"].append("expression/{sample}/{sample}.txt.fpkm".format(sample = sample))
    
    y["expression_samples"] = input_expression

    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))
    
def main(genomon_conf, run_conf, sample_conf):
    
    # preparation
    import genomon_pipeline.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq_single, sample_conf.bam_tofastq_pair)
    setup.create_directories(genomon_conf, run_conf, input_stages, 'rna/data/snakefile.txt')
    bam_tofastq_stages = (sample_conf.bam_tofastq_single, sample_conf.bam_tofastq_pair)
    setup.touch_bam_tofastq(genomon_conf, run_conf, bam_tofastq_stages)
    
    # link fastq
    linked_fastqs = setup.link_input_fastq(genomon_conf, run_conf, sample_conf.fastq, sample_conf.fastq_src)
    
    # bam import
    output_bams = setup.link_import_bam(genomon_conf, run_conf, sample_conf.bam_import, '.Aligned.sortedByCoord.out.bam', '.Aligned.sortedByCoord.out.bam.bai', "star")
    
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

    import genomon_pipeline.rna.resource.star_align as rs_star_align
    align_bams = rs_star_align.configure(genomon_conf, run_conf, sample_conf)
    output_bams.update(align_bams)

    # expression
    import genomon_pipeline.rna.resource.expression as rs_expression
    rs_expression.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    # dump conf.yaml
    dump_conf_yaml(genomon_conf, run_conf, sample_conf)
