#! /usr/bin/env python

def main(genomon_conf, run_conf, sample_conf):
    
    # preparation
    import genomon_pipeline.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq)
    setup.create_directories(genomon_conf, run_conf, input_stages, 'dna/data/snakefile.txt')
    setup.touch_bam_tofastq(genomon_conf, run_conf, (sample_conf.bam_tofastq, ))
    
    # link fastq
    linked_fastqs = setup.link_input_fastq(genomon_conf, run_conf, sample_conf.fastq, sample_conf.fastq_src)
    
    # bam import
    import genomon_pipeline.dna.resource.bwa_align as rs_align
    output_bams = setup.link_import_bam(genomon_conf, run_conf, sample_conf.bam_import, rs_align.BAM_POSTFIX, rs_align.BAI_POSTFIX)
    
    # bam to fastq
    import genomon_pipeline.dna.resource.bamtofastq as rs_bamtofastq
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
    import genomon_pipeline.dna.resource.mutation_dummy as rs_mutation
    output_mutations = rs_mutation.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    # sv
    import genomon_pipeline.dna.resource.sv_dummy as rs_sv
    output_svs = rs_sv.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    # dump conf.yaml
    y = setup.dump_yaml_input_section(
        genomon_conf, 
        run_conf, 
        (sample_conf.bam_tofastq, ),
        sample_conf.fastq,
        sample_conf.bam_import, 
        rs_align.OUTPUT_FORMAT
    )
    y["output_files"].extend(output_mutations)
    y["output_files"].extend(output_svs)
    
    y["mutation_samples"] = {}
    for (sample, control, control_panel) in sample_conf.mutation_call:
        y["mutation_samples"][sample] = rs_align.OUTPUT_FORMAT.format(sample=sample)
        
    y["sv_samples"] = {}
    for (sample, control, control_panel) in sample_conf.sv_detection:
        y["sv_samples"][sample] = rs_align.OUTPUT_FORMAT.format(sample=sample)
        
    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))
    
