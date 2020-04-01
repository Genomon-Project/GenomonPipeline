# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import os
import sys
import shutil
import unittest
import genomon_pipeline.dna.sample_conf as sc

class SubmitTest(unittest.TestCase):
    
    SAMPLE_DIR = "/tmp/temp-test/genomon_test_dna_samplesheet"
    REMOVE = True
    
    # init class
    @classmethod
    def setUpClass(self):
        os.makedirs(self.SAMPLE_DIR, exist_ok = True)
        touch_files = [
            "/A1.fastq",
            "/A2.fastq",
            "/B1.fq",
            "/B2.fq",
            "/C1_1.fq",
            "/C1_2.fq",
            "/C2_1.fq",
            "/C2_2.fq",
            "/A.markdup.bam",
            "/A.markdup.bam.bai",
            "/B.markdup.bam",
            "/B.markdup.bai",
        ]
        for p in touch_files:
            open(self.SAMPLE_DIR + p, "w").close()
    
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
    
    def test1_01_allinone(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
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
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.fastq, {
            'A_tumor': [[self.SAMPLE_DIR + '/A1.fastq'], [self.SAMPLE_DIR + '/A2.fastq']], 
            'A_tumor_2': [[self.SAMPLE_DIR + '/A1.fastq',self.SAMPLE_DIR + '/A1.fastq'], [self.SAMPLE_DIR + '/A2.fastq',self.SAMPLE_DIR + '/A2.fastq']], 
            'pool1': [[self.SAMPLE_DIR + '/B1.fq'], [self.SAMPLE_DIR + '/B2.fq']], 
            'pool2': [[self.SAMPLE_DIR + '/C1_1.fq', self.SAMPLE_DIR + '/C1_2.fq'], [self.SAMPLE_DIR + '/C2_1.fq', self.SAMPLE_DIR + '/C2_2.fq']]
        })
        #import pprint
        #pprint.pprint(sample_conf.fastq_src)
        self.assertEqual(sample_conf.fastq_src, {
            'A_tumor': [self.SAMPLE_DIR + '/A1.fastq', self.SAMPLE_DIR + '/A2.fastq'], 
            'A_tumor_2': [self.SAMPLE_DIR + '/A1.fastq', self.SAMPLE_DIR + '/A2.fastq', self.SAMPLE_DIR + '/A1.fastq', self.SAMPLE_DIR + '/A2.fastq'], 
            'pool1': [self.SAMPLE_DIR + '/B1.fq', self.SAMPLE_DIR + '/B2.fq'], 
            'pool2': [self.SAMPLE_DIR + '/C1_1.fq', self.SAMPLE_DIR + '/C2_1.fq', self.SAMPLE_DIR + '/C1_2.fq', self.SAMPLE_DIR + '/C2_2.fq']
        })
        
        #pprint.pprint(sample_conf.bam_tofastq)
        self.assertEqual(sample_conf.bam_tofastq, {
            'A_control': self.SAMPLE_DIR + '/A.markdup.bam',
            'A_control_2': self.SAMPLE_DIR + '/A.markdup.bam;' + self.SAMPLE_DIR + '/A.markdup.bam'
        })
        #pprint.pprint(sample_conf.bam_tofastq_src)
        self.assertEqual(sample_conf.bam_tofastq_src, {
            'A_control': [self.SAMPLE_DIR + '/A.markdup.bam'],
            'A_control_2': [self.SAMPLE_DIR + '/A.markdup.bam', self.SAMPLE_DIR + '/A.markdup.bam'],
        })
        self.assertEqual(sample_conf.bam_import, {'pool3': self.SAMPLE_DIR + '/B.markdup.bam'})
        self.assertEqual(sample_conf.bam_import_src, {'pool3': [self.SAMPLE_DIR + '/B.markdup.bam', self.SAMPLE_DIR + '/B.markdup.bai']})
        self.assertEqual(sample_conf.mutation_call, [('A_tumor', 'A_control', 'list1'), ('A_control', None, None)])
        self.assertEqual(sample_conf.sv_detection, [('A_tumor', 'A_control', 'list1')])
        self.assertEqual(sample_conf.qc, ['A_tumor', 'A_tumor_2', 'A_control', 'A_control_2', 'pool1', 'pool2', 'pool3'])
        self.assertEqual(sample_conf.control_panel, {'list1': ['pool1', 'pool2', 'pool3']})

    # --------------------------------------------------------------------
    # ok
    # --------------------------------------------------------------------
    def test1_02_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
A_control,{sample_dir}/A.markdup.bam
pool1,{sample_dir}/A.markdup.bam

[mutation_call]
A_tumor

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.mutation_call, [('A_tumor', None, None)])

    def test1_03_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
A_control,{sample_dir}/A.markdup.bam
pool1,{sample_dir}/A.markdup.bam

[mutation_call]
A_tumor,A_control

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.mutation_call, [('A_tumor', 'A_control', None)])

    def test1_04_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
A_control,{sample_dir}/A.markdup.bam
pool1,{sample_dir}/A.markdup.bam
pool2,{sample_dir}/A.markdup.bam

[mutation_call]
A_tumor,A_control,list1,list2

[controlpanel]
list1,pool1
list2,pool2
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.mutation_call, [('A_tumor', 'A_control', 'list1')])

    def test1_05_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
A_control,{sample_dir}/A.markdup.bam
pool1,{sample_dir}/A.markdup.bam

[sv_detection]
A_tumor

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.sv_detection, [('A_tumor', None, None)])

    def test1_06_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
A_control,{sample_dir}/A.markdup.bam
pool1,{sample_dir}/A.markdup.bam

[sv_detection]
A_tumor,A_control

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.sv_detection, [('A_tumor', 'A_control', None)])


    def test1_07_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam
A_control,{sample_dir}/A.markdup.bam
pool1,{sample_dir}/A.markdup.bam
pool2,{sample_dir}/A.markdup.bam

[sv_detection]
A_tumor,A_control,list1,list2

[controlpanel]
list1,pool1
list2,pool2
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.sv_detection, [('A_tumor', 'A_control', 'list1')])

    # --------------------------------------------------------------------
    # Not Exist
    # --------------------------------------------------------------------
    def test2_01_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A3.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_02_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B3.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_03_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_3.fq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_04_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_tofastq],,,,
A_control,{sample_dir}/X.markdup.bam,,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_05_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_import],,,,
pool3,{sample_dir}/X.markdup.bam,,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    # --------------------------------------------------------------------
    # Not defined
    # --------------------------------------------------------------------
    def test3_01_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[mutation_call]
B_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_02_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
A_control,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[mutation_call]
A_tumor,B_control,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_03_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
A_control,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[mutation_call]
A_tumor,A_control,list2

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_04_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[sv_detection]
B_tumor,None,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_05_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
A_control,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[sv_detection]
A_tumor,B_control,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_06_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
A_control,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[sv_detection]
A_tumor,A_control,list2

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_07_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[qc]
B_tumor
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_08_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
pool1,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[controlpanel]
list1,pool100
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    # --------------------------------------------------------------------
    # Duplicate
    # --------------------------------------------------------------------
    def test4_01_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_02_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_tofastq]
A_control,{sample_dir}/A.markdup.bam
A_control,{sample_dir}/A.markdup.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_03_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_import],,,,
pool3,{sample_dir}/B.markdup.bam,,,
pool3,{sample_dir}/B.markdup.bam,,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_04_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
[bam_tofastq],,,,
A_tumor,{sample_dir}/A.markdup.bam,,,
[bam_import],,,,
pool1,{sample_dir}/B.markdup.bam,,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_05_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
[bam_tofastq],,,,
pool3,{sample_dir}/A.markdup.bam,,,
[bam_import],,,,
pool3,{sample_dir}/B.markdup.bam,,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_06_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
pool3,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
[bam_tofastq],,,,
A_tumor,{sample_dir}/A.markdup.bam,,,
[bam_import],,,,
pool3,{sample_dir}/B.markdup.bam,,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)


    def test4_07_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,

[bam_tofastq],,,,
A_control,{sample_dir}/A.markdup.bam,,,

[mutation_call],,,,
A_tumor,A_control,list1,,
A_tumor,None,None,,

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_08_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,

[bam_tofastq],,,,
A_control,{sample_dir}/A.markdup.bam,,,

[sv_detection],,,,
A_tumor,A_control,list1,,
A_tumor,None,None,,

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_09_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[qc]
A_tumor
A_tumor
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_10_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
pool2,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,

[mutation_call]
A_tumor,None,list1

[controlpanel]
list1,pool1
list1,pool2
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    # --------------------------------------------------------------------
    # Format error
    # --------------------------------------------------------------------
    def test5_01_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_02_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,{sample_dir}/A2.fastq
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_03_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_tofastq]
{sample_dir}/A.markdup.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_04_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_tofastq]
A_tumor,{sample_dir}/A.markdup.bam,{sample_dir}/B.markdup.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_05_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_import]
{sample_dir}/B.markdup.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_06_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_import]
A_tumor,{sample_dir}/A.markdup.bam,{sample_dir}/B.markdup.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_07_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[mutation_call]
A_tumor,None,list1

[controlpanel]
list1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sample_conf = sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

if __name__ == '__main__':
    unittest.main()

