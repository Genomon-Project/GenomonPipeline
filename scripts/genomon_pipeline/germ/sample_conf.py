#! /usr/bin/env python
import genomon_pipeline.core.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):
    SECTION_FASTQ = "fastq"
    SECTION_BAM_IMPORT = "bam-import"
    SECTION_BAM_TOFASTQ = "bam-tofastq"
    SECTION_HTCALL = "gatk-haplotypecaller-parabrics-compatible"
    
    def __init__(self, sample_conf_file, exist_check = True):

        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq = {}
        self.bam_tofastq_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.haplotype_call = []
        self.exist_check = exist_check
        
        self.parse_file(sample_conf_file)

    def parse_data(self, _data):
        
        input_sections = [self.SECTION_FASTQ, self.SECTION_BAM_IMPORT, self.SECTION_BAM_TOFASTQ]
        analysis_sections = [self.SECTION_HTCALL]
        controlpanel_sections = []
        splited = self.split_section_data(_data, input_sections, analysis_sections, controlpanel_sections)
        
        sample_ids = []
        if self.SECTION_FASTQ in splited:
            parsed_fastq = self.parse_data_fastq_pair(splited[self.SECTION_FASTQ])
            self.fastq.update(parsed_fastq["fastq"])
            self.fastq_src.update(parsed_fastq["fastq_src"])
            sample_ids += parsed_fastq["fastq"].keys()
            
        if self.SECTION_BAM_TOFASTQ in splited:
            parsed_bam_tofastq = self.parse_data_bam_tofastq(splited[self.SECTION_BAM_TOFASTQ])
            self.bam_tofastq.update(parsed_bam_tofastq["bam_tofastq"])
            self.bam_tofastq_src.update(parsed_bam_tofastq["bam_tofastq_src"])
            sample_ids += parsed_bam_tofastq["bam_tofastq"].keys()
            
        if self.SECTION_BAM_IMPORT in splited:
            parsed_bam_import = self.parse_data_bam_import(splited[self.SECTION_BAM_IMPORT])
            self.bam_import.update(parsed_bam_import["bam_import"])
            self.bam_import_src.update(parsed_bam_import["bam_import_src"])
            sample_ids += parsed_bam_import["bam_import"].keys()
            
        if self.SECTION_HTCALL in splited:
            self.haplotype_call += self.parse_data_general(splited[self.SECTION_HTCALL])
        