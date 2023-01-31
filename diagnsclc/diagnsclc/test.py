import unittest
import input
import apply
import encoder
import numpy as np
import logging

class FileNameTest(unittest.TestCase):
    
    def test_fastq(self):
        f = 'SRR123.fastq'
        self.assertEqual(input.get_file_name(f, '.fastq'), 'SRR123')
    
    def test_bam(self):
        f = 'SRR123.bam'
        self.assertEqual(input.get_file_name(f, '.bam'), 'SRR123')
        
    def test_other(self):
        f = 'SRR123.txt'
        self.assertNotEqual(input.get_file_name(f, 'txt'), 'SRR123')
        
        
class MainTest(unittest.TestCase):
    
    def test_main_one_argv(self):
        arg = 'SRR7785623.fastq'
        try:
            input.main(['python_file.py', arg])
            self.assertTrue(True)
        except:
            self.assertTrue(False)
            
class EncoderTest(unittest.TestCase):
    
    def test_KMer_encoder(self):
        s = 'ATGC'
        self.assertIs(type(encoder.K_mer(s)), np.ndarray)

# class ApplyTest(unittest.TestCase):
    
#     def test_create_patient_file(self):
#         with self.assertLogs() as captured:
#             apply.apply()
    
if __name__ == '__main__':
    unittest.main()