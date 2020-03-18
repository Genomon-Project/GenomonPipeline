#! /usr/bin/env python
import genomon_pipeline.core.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):

    def __init__(self, sample_conf_file, no_exist_check = False):
        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq = {}
        self.bam_tofastq_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.control_panel = {}
        self.fusion = []
        self.expression = []
        self.intron_retention = []
        self.qc = []
        self.no_exist_check = no_exist_check
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
        
        input_sections = ["[fastq]", "[bam_import]", "[bam_tofastq]"]
        analysis_sections = ["[fusion]", "[expression]", "[ir_count]", "[qc]"]
        controlpanel_sections = ["[controlpanel]"]
        splited = self.split_section_data(_data, input_sections, analysis_sections, controlpanel_sections)
        
        if "[fastq]" in splited:
            parsed_fastq = self.parse_data_fastq(self, splited["[fastq]"])
            self.fastq.update(parsed_fastq["fastq"])
            self.fastq_src.update(parsed_fastq["fastq_src"])
        
        if "[bam_tofastq]" in splited:
            parsed_bam_tofastq = self.parse_data_bam_tofastq(self, splited["[bam_tofastq]"])
            self.bam_tofastq.update(parsed_bam_tofastq["bam_tofastq"])
            self.bam_tofastq_src.update(parsed_bam_tofastq["bam_tofastq_src"])
        
        if "[bam_import]" in splited:
            parsed_bam_import = self.parse_data_bam_import(self, splited["[bam_import]"])
            self.bam_import.update(parsed_bam_import["bam_import"])
            self.bam_import_src.update(parsed_bam_import["bam_import_src"])
        
        if "[expression]" in splited:
            self.expression.extend(self.parse_data_general(self, splited["[expression]"]))
        
        if "[ir_count]" in splited:
            self.intron_retention.extend(self.parse_data_general(self, splited["[ir_count]"]))
        
        if "[qc]" in splited:
            self.qc.extend(self.parse_data_general(self, splited["[qc]"]))
        
        if "[controlpanel]" in splited:
            self.control_panel.update(self.parse_data_controlpanel(self, splited["[controlpanel]"]))
        
        if "[fusion]" in splited:
            self.fusion.extend(self.parse_data_fusion(self, splited["[fusion]"], self.control_panel.keys()))
