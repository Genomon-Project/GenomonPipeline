#! /usr/bin/env python

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
    import genomon_pipeline.resource.setup_common as setup
    setup.create_directories(genomon_conf, run_conf, sample_conf, 'snakefile_dna')
    setup.link_input_fastq(genomon_conf, run_conf, sample_conf)
    output_bams = setup.link_import_bam(genomon_conf, run_conf, sample_conf, '.markdup.bam', '.markdup.bam.bai')
    
    import genomon_pipeline.resource.bwa_align
    align_bams = genomon_pipeline.resource.bwa_align.configure(genomon_conf, run_conf, sample_conf)
    output_bams.update(align_bams)
    
    import genomon_pipeline.resource.mutation_dummy
    genomon_pipeline.resource.mutation_dummy.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    import genomon_pipeline.resource.sv_dummy
    genomon_pipeline.resource.sv_dummy.configure(output_bams, genomon_conf, run_conf, sample_conf)
    
    dump_conf_yaml(genomon_conf, run_conf, sample_conf)
