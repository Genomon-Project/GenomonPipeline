# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import os
import sys
import shutil
import unittest
import genomon_pipeline.rna.sample_conf as sc

class SubmitTest(unittest.TestCase):
    
    SAMPLE_DIR = "/tmp/temp-test/genomon_test_rna_samplesheet"
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
            "/A.Aligned.sortedByCoord.out.bam",
            "/A.Aligned.sortedByCoord.out.bam.bai",
            "/B.Aligned.sortedByCoord.out.bam",
            "/B.Aligned.sortedByCoord.out.bai",
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
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.fastq, {
            'A_tumor': [[self.SAMPLE_DIR + '/A1.fastq'], [self.SAMPLE_DIR + '/A2.fastq']], 
            'A_tumor_s': [[self.SAMPLE_DIR + '/A1.fastq']], 
            'pool1': [[self.SAMPLE_DIR + '/B1.fq'], [self.SAMPLE_DIR + '/B2.fq']], 
            'pool1_s': [[self.SAMPLE_DIR + '/B1.fq']], 
            'pool2': [[self.SAMPLE_DIR + '/C1_1.fq', self.SAMPLE_DIR + '/C1_2.fq'], [self.SAMPLE_DIR + '/C2_1.fq', self.SAMPLE_DIR + '/C2_2.fq']],
            'pool2_s': [[self.SAMPLE_DIR + '/C1_1.fq', self.SAMPLE_DIR + '/C1_2.fq']],
        })
        
        self.assertEqual(sample_conf.fastq_src, {
            'A_tumor': [self.SAMPLE_DIR + '/A1.fastq', self.SAMPLE_DIR + '/A2.fastq'], 
            'A_tumor_s': [self.SAMPLE_DIR + '/A1.fastq'], 
            'pool1': [self.SAMPLE_DIR + '/B1.fq', self.SAMPLE_DIR + '/B2.fq'], 
            'pool1_s': [self.SAMPLE_DIR + '/B1.fq'], 
            'pool2': [self.SAMPLE_DIR + '/C1_1.fq', self.SAMPLE_DIR + '/C2_1.fq', self.SAMPLE_DIR + '/C1_2.fq', self.SAMPLE_DIR + '/C2_2.fq'],
            'pool2_s': [self.SAMPLE_DIR + '/C1_1.fq', self.SAMPLE_DIR + '/C1_2.fq'],
        })

        self.assertEqual(sample_conf.bam_tofastq_pair, {'A_control': self.SAMPLE_DIR + '/A.Aligned.sortedByCoord.out.bam'})
        self.assertEqual(sample_conf.bam_tofastq_pair_src, {'A_control': [self.SAMPLE_DIR + '/A.Aligned.sortedByCoord.out.bam']})
        self.assertEqual(sample_conf.bam_tofastq_single, {'A_control_s': self.SAMPLE_DIR + '/A.Aligned.sortedByCoord.out.bam'})
        self.assertEqual(sample_conf.bam_tofastq_single_src, {'A_control_s': [self.SAMPLE_DIR + '/A.Aligned.sortedByCoord.out.bam']})
        self.assertEqual(sample_conf.bam_import, {'pool3': self.SAMPLE_DIR + '/B.Aligned.sortedByCoord.out.bam'})
        self.assertEqual(sample_conf.bam_import_src, {'pool3': [self.SAMPLE_DIR + '/B.Aligned.sortedByCoord.out.bam', self.SAMPLE_DIR + '/B.Aligned.sortedByCoord.out.bai']})
        self.assertEqual(sample_conf.fusion, [('A_tumor', 'list1'), ('A_control', None)])
        self.assertEqual(sample_conf.expression, ['A_tumor', 'A_tumor_s'])
        self.assertEqual(sample_conf.qc, ['A_tumor', 'A_control', 'pool1', 'pool2', 'pool3'])
        self.assertEqual(sample_conf.control_panel, {'list1': ['pool1', 'pool2', 'pool3']})

    # --------------------------------------------------------------------
    # ok
    # --------------------------------------------------------------------
    def test1_02_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq_single]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
A_control,{sample_dir}/A.Aligned.sortedByCoord.out.bam
pool1,{sample_dir}/A.Aligned.sortedByCoord.out.bam

[fusion]
A_tumor

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.fusion, [('A_tumor', None)])

    def test1_03_ok(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """
[bam_tofastq_single]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam
A_control,{sample_dir}/A.Aligned.sortedByCoord.out.bam
pool1,{sample_dir}/A.Aligned.sortedByCoord.out.bam

[fusion]
A_tumor,list1,list2

[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.fusion, [('A_tumor', 'list1')])

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
        data = """[bam_tofastq_single],,,,
A_control,{sample_dir}/X.Aligned.sortedByCoord.out.bam,,,
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
pool3,{sample_dir}/X.Aligned.sortedByCoord.out.bam,,,
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

[fusion]
B_tumor,None
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
pool1,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,

[fusion]
A_tumor,list2

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

    def test3_03_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[expression]
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

    def test3_04_undefine(self):
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

    def test3_05_undefine(self):
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
        data = """[bam_tofastq_single]
A_control,{sample_dir}/A.Aligned.sortedByCoord.out.bam
A_control,{sample_dir}/A.Aligned.sortedByCoord.out.bam
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
pool3,{sample_dir}/B.Aligned.sortedByCoord.out.bam,,,
pool3,{sample_dir}/B.Aligned.sortedByCoord.out.bam,,,
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
[bam_tofastq_single],,,,
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam,,,
[bam_import],,,,
pool1,{sample_dir}/B.Aligned.sortedByCoord.out.bam,,,
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
[bam_tofastq_single],,,,
pool3,{sample_dir}/A.Aligned.sortedByCoord.out.bam,,,
[bam_import],,,,
pool3,{sample_dir}/B.Aligned.sortedByCoord.out.bam,,,
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
[bam_tofastq_single],,,,
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam,,,
[bam_import],,,,
pool3,{sample_dir}/B.Aligned.sortedByCoord.out.bam,,,
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

[bam_tofastq_single],,,,
A_control,{sample_dir}/A.Aligned.sortedByCoord.out.bam,,,

[fusion],,,,
A_tumor,list1,,
A_tumor,None,,

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
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[expression]
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

[fusion]
A_tumor,list1

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
    def test5_01_unformat_ok(self):
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

        self.assertFalse(fail)

    def test5_02_unformat_ok(self):
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

        self.assertFalse(fail)

    def test5_03_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[bam_tofastq_single]
{sample_dir}/A.Aligned.sortedByCoord.out.bam
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
        data = """[bam_tofastq_single]
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam,{sample_dir}/B.Aligned.sortedByCoord.out.bam
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
{sample_dir}/B.Aligned.sortedByCoord.out.bam
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
A_tumor,{sample_dir}/A.Aligned.sortedByCoord.out.bam,{sample_dir}/B.Aligned.sortedByCoord.out.bam
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

[fusion]
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

