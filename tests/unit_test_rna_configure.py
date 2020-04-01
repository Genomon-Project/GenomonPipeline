# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import os
import sys
import shutil
import unittest
import subprocess
import snakemake

def func_path (root, name):
    wdir = root + "/" + name
    ss_path = root + "/" + name + ".csv"
    return (wdir, ss_path)

class ConfigureTest(unittest.TestCase):
    
    DATA_DIR = "/tmp/temp-test/genomon_test_rna_configure"
    SAMPLE_DIR = DATA_DIR + "/samples"
    REMOVE = False
    SS_NAME = "/test.csv"
    GC_NAME = "/genomon.cfg"

    # init class
    @classmethod
    def setUpClass(self):
        os.makedirs(self.SAMPLE_DIR, exist_ok = True)
        os.makedirs(self.DATA_DIR + "/reference", exist_ok = True)
        os.makedirs(self.DATA_DIR + "/image", exist_ok = True)
        touch_files = [
            "/samples/A1.fastq",
            "/samples/A2.fastq",
            "/samples/B1.fq",
            "/samples/B2.fq",
            "/samples/C1_1.fq",
            "/samples/C1_2.fq",
            "/samples/C2_1.fq",
            "/samples/C2_2.fq",
            "/samples/A.Aligned.sortedByCoord.out.bam",
            "/samples/A.Aligned.sortedByCoord.out.bam.bai",
            "/samples/B.Aligned.sortedByCoord.out.bam",
            "/samples/B.Aligned.sortedByCoord.out.bai",
            "/reference/XXX.fa",
            "/image/YYY.simg",
            "/reference/ZZZ.gtf"
        ]
        for p in touch_files:
            open(self.DATA_DIR + p, "w").close()
        
        data_sample = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq,,
A_tumor_s,{sample_dir}/A1.fastq
pool1_s,{sample_dir}/B1.fq
pool2_s,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,,,,

[bam_tofastq_pair],,,,
A_control,{sample_dir}/A.Aligned.sortedByCoord.out.bam,,,

[bam_tofastq_single],,,,
A_control_s,{sample_dir}/A.Aligned.sortedByCoord.out.bam,,,

[bam_import],,,,
pool3,{sample_dir}/B.Aligned.sortedByCoord.out.bam,,,
,,,,
[fusion],,,,
A_tumor,list1,,
A_control,None,,
,,,,
[expression],,,,
A_tumor,,,,
A_tumor_s,,,,
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
        
        f = open(self.DATA_DIR + self.SS_NAME, "w")
        f.write(data_sample)
        f.close()

        data_conf = """[bam_tofastq]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 
params = collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY tryoq=0

[star_alignment]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 
star_option = --runThreadN 6 --outSAMstrandField intronMotif --outSAMunmapped Within --alignMatesGapMax 500000 --alignIntronMax 500000 --alignSJstitchMismatchNmax -1 -1 -1 -1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --outSAMtype BAM Unsorted
star_genome = {sample_dir}/reference
remove_fastq = True

[expression]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 
gtf = {sample_dir}/reference/ZZZ.gtf
""".format(sample_dir = self.DATA_DIR)
        
        f = open(self.DATA_DIR + self.GC_NAME, "w")
        f.write(data_conf)
        f.close()

    # terminated class
    @classmethod
    def tearDownClass(self):
        if self.REMOVE:
            shutil.rmtree(self.DATA_DIR)

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
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_02_configure(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
            "--disable_drmaa",
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_01_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_02_1_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_02_2_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_pair]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_03_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        

        data_sample = """[bam_import]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test4_01_fusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[fusion]
A_tumor,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_02_fusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
[fusion]
A_tumor,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_03_fusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
[fusion]
A_tumor,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test5_01_exp_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[expression]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_02_exp_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
[expression]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_03_exp_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
[expression]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


if __name__ == '__main__':
    unittest.main()
