#! /usr/bin/env python

def dump_conf_yaml(genomon_conf, run_conf, sample_conf):
    samples = []
    outputs =[]
    for sample in sample_conf.bam_tofastq:
        samples.append(sample)
        outputs.append("fastq/{sample}/1_1.fastq".format(sample = sample))
        outputs.append("fastq/{sample}/1_2.fastq".format(sample = sample))

    input_bwa = {}
    for sample in sample_conf.fastq:
        input_bwa[sample] = "fastq/%s/%s" % (sample, sample_conf.fastq[sample][0][0].split("/")[-1])
        input_bwa[sample] = "fastq/%s/%s" % (sample, sample_conf.fastq[sample][1][0].split("/")[-1])
        outputs.append("bam/{sample}/{sample}.markdup.bam".format(sample = sample))
    
    input_mutation = {}
    for (sample, control, control_panel) in sample_conf.mutation_call:
        input_mutation[sample] = "bam/{sample}/{sample}.markdup.bam".format(sample = sample)
        outputs.append("mutation/{sample}/{sample}.txt".format(sample = sample))
    
    input_sv = {}
    for (sample, control, control_panel) in sample_conf.sv_detection:
        input_sv[sample] = "bam/{sample}/{sample}.markdup.bam".format(sample = sample)
        outputs.append("sv/{sample}/{sample}.txt".format(sample = sample))
    
    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump({
        "samples": samples,
        "bwa_samples": input_bwa,
        "mutation_samples": input_mutation,
        "sv_samples": input_sv,
        "output_files": outputs
    }))
    
def main(genomon_conf, run_conf, sample_conf):
    
    # preparation
    import genomon_pipeline.core.setup_common as setup
    setup.create_directories(genomon_conf, run_conf, sample_conf, 'dna/data/snakefile.txt')
    setup.touch_bam_tofastq(genomon_conf, run_conf, sample_conf)
    
    # link fastq
    linked_fastqs = setup.link_input_fastq(genomon_conf, run_conf, sample_conf)
    
    # bam import
    output_bams = setup.link_import_bam(genomon_conf, run_conf, sample_conf, '.markdup.bam', '.markdup.bam.bai')
    
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
