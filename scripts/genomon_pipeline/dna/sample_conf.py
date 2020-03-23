#! /usr/bin/env python
import genomon_pipeline.core.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):
    SECTION_FASTQ = "fastq"
    SECTION_BAM_IMPORT = "bam_import"
    SECTION_BAM_TOFASTQ = "bam_tofastq"
    SECTION_MUTATION = "mutation_call"
    SECTION_SV = "sv_detection"
    SECTION_QC = "qc"
    SECTION_CONTROL_PANEL = "controlpanel"
    
    def __init__(self, sample_conf_file, exist_check = True):

        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq = {}
        self.bam_tofastq_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.mutation_call = []
        self.sv_detection = []
        self.qc = []
        self.control_panel = {}
        self.exist_check = exist_check
        
        self.parse_file(sample_conf_file)

    def parse_data_analysis(self, _data, mode, sampleID_list, controlpanel_list):
        
        analysis = []
        for row in _data:
    
            tumorID = row[0]
            if tumorID not in sampleID_list:
                err_msg = "%s section, %s is not defined" % (mode, tumorID)
                raise ValueError(err_msg)

            normalID = row[1] if len(row) >= 2 and row[1] not in ['', 'None'] else None
            if normalID is not None and normalID not in sampleID_list:
                err_msg = "%s section, %s is not defined" % (mode, normalID)
                raise ValueError(err_msg)

            controlpanelID = row[2] if len(row) >= 3 and row[2] not in ['', 'None'] else None
            if controlpanelID is not None and controlpanelID not in controlpanel_list:
                err_msg = "%s section, %s is not defined" % (mode, controlpanelID)
                raise ValueError(err_msg)

            analysis.append((tumorID, normalID, controlpanelID))
    
        return analysis
    
    def parse_data(self, _data):
        
        input_sections = [self.SECTION_FASTQ, self.SECTION_BAM_IMPORT, self.SECTION_BAM_TOFASTQ]
        analysis_sections = [self.SECTION_MUTATION, self.SECTION_SV, self.SECTION_QC]
        controlpanel_sections = [self.SECTION_CONTROL_PANEL]
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
            
        if self.SECTION_QC in splited:
            self.qc += self.parse_data_general(splited[self.SECTION_QC])
        
        if self.SECTION_CONTROL_PANEL in splited:
            self.control_panel.update(self.parse_data_controlpanel(splited[self.SECTION_CONTROL_PANEL]))
        
        if self.SECTION_MUTATION in splited:
            self.mutation_call += self.parse_data_analysis(splited[self.SECTION_MUTATION], self.SECTION_MUTATION, sample_ids, self.control_panel.keys())
        
        if self.SECTION_SV in splited:
            self.sv_detection += self.parse_data_analysis(splited[self.SECTION_SV], self.SECTION_SV, sample_ids, self.control_panel.keys())
        