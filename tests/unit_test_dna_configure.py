# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import os
import shutil
import unittest
import subprocess

class SampleSheetTest(unittest.TestCase):
    
    SAMPLE_DIR = "/tmp/genomon_test_dna_configure/"
    REMOVE = False
    SS_NAME = "test.csv"
    GC_NAME = "genomon.cfg"

    # init class
    @classmethod
    def setUpClass(self):
        os.makedirs(self.SAMPLE_DIR, exist_ok = True)
        touch_files = [
            "A1.fastq",
            "A2.fastq",
            "B1.fq",
            "B2.fq",
            "C1_1.fq",
            "C1_2.fq",
            "C2_1.fq",
            "C2_2.fq",
            "A.markdup.bam",
            "A.markdup.bam.bai",
            "B.markdup.bam",
            "B.markdup.bai",
            "XXX.fa",
            "YYY.simg"
        ]
        for p in touch_files:
            open(self.SAMPLE_DIR + p, "w").close()
        
        data_sample = """[fastq],,,,
A_tumor,{sample_dir}A1.fastq,{sample_dir}A2.fastq,,
pool1,{sample_dir}B1.fq,{sample_dir}B2.fq,,
pool2,{sample_dir}C1_1.fq;{sample_dir}C1_2.fq,{sample_dir}C2_1.fq;{sample_dir}C2_2.fq,,
,,,,
[bam_tofastq],,,,
A_control,{sample_dir}A.markdup.bam,,,
[bam_import],,,,
pool3,{sample_dir}B.markdup.bam,,,
,,,,
[mutation_call],,,,
A_tumor,A_control,list1,,
A_control,None,None,,
,,,,
[sv_detection],,,,
A_tumor,A_control,list1,,
,,,,
[qc],,,,
A_tumor,,,,
A_control,,,,
pool1,,,,
pool2,,,,
pool3,,,,
,,,,
[controlpanel]
list1,pool1,pool2,pool3
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(self.SAMPLE_DIR + self.SS_NAME, "w")
        f.write(data_sample)
        f.close()

        data_conf = """[bam_tofastq]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}YYY.simg
singularity_option = 

params = collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY tryoq=0

[bwa_alignment]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}YYY.simg
singularity_option = 

bamtofastq_option = collate=1 combs=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY tryoq=1
bwa_option = -t 8 -T 0
bwa_reference = {sample_dir}XXX.fa
bamsort_option = index=1 level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindentonly=1
bammarkduplicates_option = markthreads=2 rewritebam=1 rewritebamlevel=1 index=1 md5=1

[mutation_dummy]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}YYY.simg
singularity_option = 

[sv_dummy]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}YYY.simg
singularity_option = 
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(self.SAMPLE_DIR + self.GC_NAME, "w")
        f.write(data_conf)
        f.close()

    # terminated class
    @classmethod
    def tearDownClass(self):
        if self.REMOVE:
            shutil.rmtree(self.SAMPLE_DIR)

    # init method
    def setUp(self):
        pass

    # terminated method
    def tearDown(self):
        pass
    
    def test1_01_version(self):
        subprocess.check_call(['python', 'genomon_pipeline', '--version'])
    
    def test1_02_version(self):
        subprocess.check_call(['python', 'genomon_runner', '--version'])
    
    def test2_01_configure(self):
        options = [
            "dna",
            self.SAMPLE_DIR + self.SS_NAME,
            self.SAMPLE_DIR + "drmaa",
            self.SAMPLE_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
    
    def test2_02_configure(self):
        options = [
            "dna",
            self.SAMPLE_DIR + self.SS_NAME,
            self.SAMPLE_DIR + "qsub",
            self.SAMPLE_DIR + self.GC_NAME,
            "--disable_drmaa",
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)

    def test3_01_dryrun(self):
        import snakemake
        snakemake.snakemake(self.SAMPLE_DIR + 'drmaa/snakefile', workdir = self.SAMPLE_DIR + 'drmaa', dryrun = True)

if __name__ == '__main__':
    unittest.main()
