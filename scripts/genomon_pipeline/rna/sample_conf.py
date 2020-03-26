#! /usr/bin/env python
import genomon_pipeline.core.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):
    SECTION_FASTQ = "fastq"
    SECTION_BAM_IMPORT = "bam_import"
    SECTION_BAM_TOFASTQ_PAIR = "bam_tofastq_pair"
    SECTION_BAM_TOFASTQ_SINGLE = "bam_tofastq_single"
    SECTION_FUSION = "fusion"
    SECTION_EXPRESSION = "expression"
    SECTION_IR = "intron_retention"
    SECTION_QC = "qc"
    SECTION_CONTROL_PANEL = "controlpanel"
    
    def __init__(self, sample_conf_file, exist_check = True):
        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq_pair = {}
        self.bam_tofastq_pair_src = {}
        self.bam_tofastq_single = {}
        self.bam_tofastq_single_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.control_panel = {}
        self.fusion = []
        self.expression = []
        self.intron_retention = []
        self.qc = []
        self.exist_check = exist_check
        self.parse_file(sample_conf_file)
    
    def parse_data_fusion(self, _data, controlpanel_list):
        
        analysis = []
        for row in _data:
            controlpanelID = row[1] if len(row) >= 2 and row[1] not in ['', 'None'] else None
            if controlpanelID != None and not controlpanelID in controlpanel_list:
                err_msg = "[fusion] section, %s is not defined" % (controlpanelID)
                raise ValueError(err_msg)
            analysis.append((row[0], controlpanelID))
        return analysis
    
    def parse_data(self, _data):
        
        input_sections = [self.SECTION_FASTQ , self.SECTION_BAM_IMPORT, self.SECTION_BAM_TOFASTQ_PAIR, self.SECTION_BAM_TOFASTQ_SINGLE]
        analysis_sections = [self.SECTION_FUSION, self.SECTION_EXPRESSION, self.SECTION_IR, self.SECTION_QC]
        controlpanel_sections = [self.SECTION_CONTROL_PANEL]
        splited = self.split_section_data(_data, input_sections, analysis_sections, controlpanel_sections)
        
        if self.SECTION_FASTQ in splited:
            parsed_fastq = self.parse_data_fastq_mixed(splited[self.SECTION_FASTQ])
            self.fastq.update(parsed_fastq["fastq"])
            self.fastq_src.update(parsed_fastq["fastq_src"])

        if self.SECTION_BAM_TOFASTQ_PAIR in splited:
            parsed_bam_tofastq_pair = self.parse_data_bam_tofastq(splited[self.SECTION_BAM_TOFASTQ_PAIR])
            self.bam_tofastq_pair.update(parsed_bam_tofastq_pair["bam_tofastq"])
            self.bam_tofastq_pair_src.update(parsed_bam_tofastq_pair["bam_tofastq_src"])
        
        if self.SECTION_BAM_TOFASTQ_SINGLE in splited:
            parsed_bam_tofastq_single = self.parse_data_bam_tofastq(splited[self.SECTION_BAM_TOFASTQ_SINGLE])
            self.bam_tofastq_single.update(parsed_bam_tofastq_single["bam_tofastq"])
            self.bam_tofastq_single_src.update(parsed_bam_tofastq_single["bam_tofastq_src"])
        
        if self.SECTION_BAM_IMPORT in splited:
            parsed_bam_import = self.parse_data_bam_import(splited[self.SECTION_BAM_IMPORT])
            self.bam_import.update(parsed_bam_import["bam_import"])
            self.bam_import_src.update(parsed_bam_import["bam_import_src"])
        
        if self.SECTION_EXPRESSION in splited:
            self.expression += self.parse_data_general(splited[self.SECTION_EXPRESSION])
        
        if self.SECTION_IR in splited:
            self.intron_retention += self.parse_data_general(splited[self.SECTION_IR])
        
        if self.SECTION_QC in splited:
            self.qc += self.parse_data_general(splited[self.SECTION_QC])
        
        if self.SECTION_CONTROL_PANEL in splited:
            self.control_panel.update(self.parse_data_controlpanel(splited[self.SECTION_CONTROL_PANEL]))
        
        if self.SECTION_FUSION in splited:
            self.fusion += self.parse_data_fusion(splited[self.SECTION_FUSION], self.control_panel.keys())
