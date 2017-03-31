import sys
import unittest

sys.path.insert(0, "..")

import util

class UtilTests(unittest.TestCase):
    
    def test_pos2idx(self):
        pos = [0,0,10]
        image_size = [128,128,128]
        idx = util.pos2idx(pos, image_size)
        self.assertEqual(idx, 10)

    def test_idx2pos(self):
        idx = 10
        image_size = [128,128,128]
        pos = util.idx2pos(idx, image_size)
        self.assertTrue(len(pos)==3)
        self.assertEqual(pos[0], 0)
        self.assertEqual(pos[1], 0)
        self.assertEqual(pos[2], 10)



if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(UtilTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
    