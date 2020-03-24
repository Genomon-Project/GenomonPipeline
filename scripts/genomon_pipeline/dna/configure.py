#! /usr/bin/env python

def dump_conf_yaml(genomon_conf, run_conf, sample_conf):
    
    import genomon_pipeline.core.setup_common as setup
    y = setup.dump_yaml_input_section(
        genomon_conf, 
        run_conf, 
        (sample_conf.bam_tofastq),
        sample_conf.fastq,
        sample_conf.bam_import, 
        "bam/{sample}/{sample}.markdup.bam"
    )
    
    input_mutation = {}
    for (sample, control, control_panel) in sample_conf.mutation_call:
        input_mutation[sample] = "bam/{sample}/{sample}.markdup.bam".format(sample = sample)
        y["output_files"].append("mutation/{sample}/{sample}.txt".format(sample = sample))
    
    input_sv = {}
    for (sample, control, control_panel) in sample_conf.sv_detection:
        input_sv[sample] = "bam/{sample}/{sample}.markdup.bam".format(sample = sample)
        y["output_files"].append("sv/{sample}/{sample}.txt".format(sample = sample))
    
    y["mutation_samples"] = input_mutation
    y["sv_samples"] = input_sv

    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))
    
def main(genomon_conf, run_conf, sample_conf):
    
    # preparation
    import genomon_pipeline.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq)
    setup.create_directories(genomon_conf, run_conf, input_stages, 'dna/data/snakefile.txt')
    setup.touch_bam_tofastq(genomon_conf, run_conf, (sample_conf.bam_tofastq))
    
    # link fastq
    linked_fastqs = setup.link_input_fastq(genomon_conf, run_conf, sample_conf.fastq, sample_conf.fastq_src)
    
    # bam import
    output_bams = setup.link_import_bam(genomon_conf, run_conf, sample_conf.bam_import, '.markdup.bam', '.markdup.bam.bai')
    
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

    import genomon_pipeline.dna.resource.bwa_align as rs_bwa_align
    align_bams = rs_bwa_align.configure(genomon_conf, run_conf, sample_conf)
    output_bams.update(align_bams)

    # mutation
    import genomon_pipeline.dna.resource.mutation_dummy as rs_mutation
    rs_mutation.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    # sv
    import genomon_pipeline.dna.resource.sv_dummy as rs_sv
    rs_sv.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    # dump conf.yaml
    dump_conf_yaml(genomon_conf, run_conf, sample_conf)
