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
    
    DATA_DIR = "/tmp/temp-test/genomon_test_dna_configure"
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
            "/samples/A.markdup.bam",
            "/samples/A.markdup.bam.bai",
            "/samples/B.markdup.bam",
            "/samples/B.markdup.bai",
            "/reference/XXX.fa",
            "/reference/gap.txt",
            "/reference/human.hg19.genome",
            "/reference/refGene.coding.exon.151207.bed",
            "/image/YYY.simg",
        ]
        for p in touch_files:
            open(self.DATA_DIR + p, "w").close()
        
        data_sample = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
A_tumor_2,{sample_dir}/A1.fastq;{sample_dir}/A1.fastq,{sample_dir}/A2.fastq;{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq,,
,,,,
[bam_tofastq],,,,
A_control,{sample_dir}/A.markdup.bam,,,
A_control_2,{sample_dir}/A.markdup.bam;{sample_dir}/A.markdup.bam,,,

[bam_import],,,,
pool3,{sample_dir}/B.markdup.bam,,,
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
A_tumor_2,,,,
A_control,,,,
A_control_2,,,,
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

[bwa_alignment]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

bamtofastq_option = collate=1 combs=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY tryoq=1
bwa_option = -t 8 -T 0
bwa_reference = {sample_dir}/reference/XXX.fa
bamsort_option = index=1 level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindentonly=1
bammarkduplicates_option = markthreads=2 rewritebam=1 rewritebamlevel=1 index=1 md5=1
remove_fastq = True

[mutation_dummy]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

[sv_dummy]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

[qc_bamstats]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

[qc_coverage]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

coverage    = 2,10,20,30,40,50,100
wgs_flag = False
wgs_incl_bed_width = 1000000
wgs_i_bed_lines = 10000
wgs_i_bed_width = 100
samtools_params = -F 3332 -f 2
grc_flag = True
genome_size = {sample_dir}/reference/human.hg19.genome
gaptxt = {sample_dir}/reference/gap.txt
bait_file = {sample_dir}/reference/refGene.coding.exon.151207.bed

[qc_merge]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

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
            "dna",
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
            "dna",
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
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_02_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
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
A_tumor,{sample_dir}/A.markdup.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test4_01_mutation_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[mutation_call]
A_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_02_mutation_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
[mutation_call]
A_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_03_mutation_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
A_tumor,{sample_dir}/A.markdup.bam
[mutation_call]
A_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test5_01_sv_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[sv_detection]
A_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_02_sv_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
[sv_detection]
A_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_03_sv_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
A_tumor,{sample_dir}/A.markdup.bam
[sv_detection]
A_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "dna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'genomon_pipeline'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

#    def test_dag(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[bam_import]
#A_tumor,{sample_dir}/A.markdup.bam
#[sv_detection]
#A_tumor,None,None
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "dna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC_NAME,
#        ]
#        subprocess.check_call(['python', 'genomon_pipeline'] + options)
#        sys.stdout = open(wdir + "/dag.txt", "w")
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, forceall = True, printdag = True)
#        sys.stdout.close()
#        sys.stdout = sys.__stdout__
#        subprocess.check_call(('cat ' + wdir + '/dag.txt | dot -Tpng > ' + wdir + '/dag.png').split(" "))
#
#        self.assertTrue(success)

if __name__ == '__main__':
    unittest.main()
