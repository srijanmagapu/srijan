import unittest
import numpy.testing as ntest

import math

#Author - Srijan Magapu - Symphonium Project


from .. import wrapper as fingerprinter

filename = ""

# The wrapper class
class calculate_frequency_bands_test(unittest.TestCase):
  
  # Test the length of the frequency band boundary array. Should 
  # be the number of bits/fingerprint + 2
  def test_different_bitsizes(self):
    f = fingerprinter.Fingerprinter(filepath=filename, framewidth=0.37, overlap=0.1)
    
    for x in range(1, 33):
      f.FINGERPRINT_NBITS = x
      self.assertEqual((f.FINGERPRINT_NBITS + 2,), f.calculate_frequency_bands(4097, 0.37).shape)
    
    f.close()


  
  
  # Just some random value compared with what the function outputs
  def test_default(self):
    f = fingerprinter.Fingerprinter(filepath=filename, framewidth=0.37, overlap=0.1)
    try:
      ntest.assert_array_equal(f.calculate_frequency_bands(4079, 0.37), [111,118,125,132,140,148,157,166,176,186,197,209,221,234,248,263,278,295,312,331,350,371,393,416,441,467,495,524,555,588,623,660,699,740])
    except AssertionError:
      self.fail("Frequency bands are not being calculated correctly.")


if __name__ == "__main__":
  unittest.main()
