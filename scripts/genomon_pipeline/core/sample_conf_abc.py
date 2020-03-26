#! /usr/bin/env python
import os

class Sample_conf_abc(object):

    def __init__(self, sample_conf_file, exist_check = True):
        pass
        
    def parse_file(self, file_path):

        file_ext = file_path.split(".")[-1].lower()

        file_data = []
        if file_ext == 'csv':
            file_data = self.parse_csv(file_path)
        elif file_ext in ['txt', 'tsv']:
            file_data = self.parse_tsv(file_path)
        else:
            raise NotImplementedError("currently, we can just accept tsv and csv formats")
 

        file_data_trimmed = []
        for line_data in file_data:
       
            # skip empty lines
            if len(line_data) == 0: continue
 
            # line starting with '#' is comment
            if line_data[0].startswith('#'): continue
             
            # remove spaces
            line_data = list(map(lambda x: x.strip(' '), line_data))

            # skip if all the elements are empty
            if len(line_data) == line_data.count(''): continue

            file_data_trimmed.append(line_data)


        self.parse_data(file_data_trimmed)


    def parse_csv(self, file_path):

        _file_data = []
        import csv
        with open(file_path, 'r') as hIN:
            csv_obj = csv.reader(hIN)
            for cells in csv_obj:
                tempdata = []
                row_len = 0
                for cell in cells:
                    row_len += len(cell)
                    if (len(cell) == 0) and (row_len > 0):
                        continue
                    tempdata.append(cell)
                
                if row_len > 0:
                    _file_data.append(tempdata)

        return _file_data


    def parse_tsv(self, file_path):

        _file_data = []
        with open(file_path, 'r') as hIN:
            for line in hIN:
                F = line.rstrip().split('\t')
                tempdata = []
                row_len = 0
                for cell in F:
                    row_len += len(cell)
                    if (len(cell) == 0) and (row_len > 0):
                        continue
                    tempdata.append(cell)
                
                if row_len > 0:
                    _file_data.append(tempdata)

        return _file_data

    def _link_sources(self, file_path):
        
        links = []
        link_dst = os.path.abspath(file_path)
        while(1):
            if os.path.islink(link_dst):
                link_src = os.readlink(link_dst)
                if link_src.startswith("/"):
                    links.append(os.path.abspath(link_src))
                else:
                    links.append(os.path.abspath(os.path.dirname(link_dst) + "/" + link_src))
                link_dst = os.path.abspath(link_src)
            else:
                break

        return links
    
    def _exists(self, file_path):
        if self.exist_check:
            return os.path.exists(file_path)
        return True
    
    def split_section_data(self, _data, input_sections, analysis_sections, controlpanel_sections = []):
        
        split_data = {}
        sampleID_list = []
        except_sections = input_sections + analysis_sections + controlpanel_sections
        mode = ""
        for row in _data:
            if row[0].startswith('['):
                mode = row[0].lower().replace("[", "").replace("]", "")
                if not mode in except_sections:
                    err_msg = "Unexcepted section %s. " % (mode)
                    err_msg += "Section name should be either of %s. "  % (",".join(except_sections))
                    err_msg += "Also, sample name should not start with '['."
                    raise ValueError(err_msg)
                    
                split_data[mode] = []
                continue
            
            sampleID = row[0]
            if sampleID == 'None':
                err_msg = "None can not be used as sampleID"
                raise ValueError(err_msg)

            if mode in input_sections:
                if sampleID in sampleID_list:
                    err_msg = sampleID + " is duplicated"
                    raise ValueError(err_msg)

                sampleID_list.append(sampleID)
            elif mode in analysis_sections:
                if not sampleID in sampleID_list:
                    err_msg = "%s section, %s is not defined" % (mode, sampleID)
                    raise ValueError(err_msg)
            
            elif mode in controlpanel_sections:
                for i, r in enumerate(row):
                    if i == 0:
                        continue
                    if not r in sampleID_list:
                        err_msg = "%s section, %s is not defined" % (mode, r)
                        raise ValueError(err_msg)
                
            split_data[mode].append(row)
        
        for mode in analysis_sections:
            if not mode in split_data:
                continue
            l = [x[0] for x in split_data[mode]]
            dup = [x for x in set(l) if l.count(x) > 1]
            if len(dup) > 0:
                err_msg = "%s section, %s is duplicated" % (mode, ",".join(dup))
                raise ValueError(err_msg)
        
        for mode in controlpanel_sections:
            if not mode in split_data:
                continue
            l = [x[0] for x in split_data[mode]]
            dup = [x for x in set(l) if l.count(x) > 1]
            if len(dup) > 0:
                err_msg = "%s section, %s is duplicated" % (mode, ",".join(dup))
                raise ValueError(err_msg)
        
        return split_data
    
    def parse_data_fastq_pair(self, _data):
        """
        pair-read
        """
        fastq = {}
        fastq_src = {}
        
        for row in _data:
            sampleID = row[0]
            
            if len(row) != 3:
                err_msg = sampleID + ": the path for read1 (and read2) should be provided"
                raise ValueError(err_msg)

            sequence1 = row[1].split(';')
            sequence2 = row[2].split(';')
            
            fastq_src[sampleID] = []
            for s in range(len(sequence1)):
                if not self._exists(sequence1[s]):
                    err_msg = sampleID + ": " + sequence1[s] +  " does not exists" 
                    raise ValueError(err_msg)
                if not self._exists(sequence2[s]):
                    err_msg = sampleID + ": " + sequence2[s] +  " does not exists" 
                    raise ValueError(err_msg)
                if sequence1[s] == sequence2[s]:
                    err_msg = sampleID + ": read1 and read2 are same path" 
                    raise ValueError(err_msg)
                
                fastq_src[sampleID].append(sequence1[s])
                fastq_src[sampleID].append(sequence2[s])
                fastq_src[sampleID].extend(self._link_sources(sequence1[s]))
                fastq_src[sampleID].extend(self._link_sources(sequence2[s]))

            fastq[sampleID] = [sequence1, sequence2]
    
        return {"fastq": fastq, "fastq_src": fastq_src}
    

    def parse_data_fastq_mixed(self, _data):
        """
        single-read or pair-read
        """
        fastq = {}
        fastq_src = {}
        
        for row in _data:
            sampleID = row[0]
            
            if len(row) < 2:
                err_msg = sampleID + ": the path for read1 (and read2) should be provided"
                raise ValueError(err_msg)
            
            pair = len(row) == 3
            sequence1 = row[1].split(';')
            if pair:
                sequence2 = row[2].split(';')
                
            fastq_src[sampleID] = []
            for s in range(len(sequence1)):
                if not self._exists(sequence1[s]):
                    err_msg = sampleID + ": " + sequence1[s] +  " does not exists" 
                    raise ValueError(err_msg)
                fastq_src[sampleID].append(sequence1[s])
                fastq_src[sampleID].extend(self._link_sources(sequence1[s]))
                                
                if pair:
                    if not self._exists(sequence2[s]):
                        err_msg = sampleID + ": " + sequence2[s] +  " does not exists" 
                        raise ValueError(err_msg)
                    fastq_src[sampleID].append(sequence2[s])
                    fastq_src[sampleID].extend(self._link_sources(sequence2[s]))
            
            if pair:
                fastq[sampleID] = [sequence1, sequence2]
            else:
                fastq[sampleID] = [sequence1]
        return {"fastq": fastq, "fastq_src": fastq_src}
    

    def parse_data_bam_tofastq(self, _data):
    
        bam_tofastq = {}
        bam_tofastq_src = {}
        
        for row in _data:
            sampleID = row[0]
            
            if len(row) != 2:
                err_msg = sampleID + ": only one bam file is allowed"
                raise ValueError(err_msg)

            bam_tofastq_src[sampleID] = []
            sequences = row[1]
            for seq in sequences.split(";"):
                if not self._exists(seq):
                    err_msg = sampleID + ": " + seq +  " does not exists"
                    raise ValueError(err_msg)
                bam_tofastq_src[sampleID].append(seq)
                bam_tofastq_src[sampleID].extend(self._link_sources(seq))

            bam_tofastq[sampleID] = sequences
        return {"bam_tofastq": bam_tofastq, "bam_tofastq_src": bam_tofastq_src}
 
    def parse_data_bam_import(self, _data):
    
        bam_import = {}
        bam_import_src = {}
        
        for row in _data:
            sampleID = row[0]
            
            if len(row) != 2:
                err_msg = sampleID + ": only one bam file is allowed"
                raise ValueError(err_msg)

            sequence = row[1]
            if not self._exists(sequence):
                err_msg = sampleID + ": " + sequence +  " does not exists"
                raise ValueError(err_msg)
            
            sequence_prefix, sequence_ext = os.path.splitext(sequence)
            if sequence_ext == ".bam":
                bam_index = ".bai"
            elif sequence_ext == ".cram":
                bam_index = ".crai"
            else:
                err_msg = sampleID + ": " + sequence +  " is unsupported file"
                raise ValueError(err_msg)

            sequence_index = ""
            if self._exists(sequence + bam_index):
                sequence_index = sequence + bam_index
            elif self._exists(sequence_prefix + bam_index):
                sequence_index = sequence_prefix + bam_index
            else:
                err_msg = sampleID + ": " + sequence +  " index does not exists"
                raise ValueError(err_msg)

            bam_import_src[sampleID] = [sequence, sequence_index]
            bam_import_src[sampleID].extend(self._link_sources(sequence))
            bam_import_src[sampleID].extend(self._link_sources(sequence_index))
            
            bam_import[sampleID] = sequence
        return {"bam_import": bam_import, "bam_import_src": bam_import_src}
    
    def parse_data_general(self, _data):
        
        sample_list = []
        for row in _data:
            sample_list.append(row[0])
        return sample_list

    def parse_data_controlpanel(self, _data):
        
        control_panel = {}
        for row in _data:
            if len(row) <= 1:
                err_msg = "[controlpanel] section, list item is none for the row: " + ','.join(row)
                raise ValueError(err_msg)

            controlpanelID = row[0]
            control_panel[controlpanelID] = row[1:]
        return control_panel
    
    
    def parse_data(self, _data ):
        pass

