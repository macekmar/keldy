#!/usr/bin/env python

from math import sqrt
import unittest

from keldy import binner

class test_binner(unittest.TestCase):

    def test_basics(self):
        b = binner.Binner_1_1_dcomplex([(-5., 5., 10)], [3])

        data = b.get_data()
        nr_values_added = b.get_nr_values_added()
        nr_values_dropped = b.get_nr_values_dropped()
        bin_times = b.get_bin_coord(0)
        bin_size = b.get_bin_size(0)

        self.assertEqual(data.shape, (10, 3))
        self.assertEqual(nr_values_dropped, 0)

class test_sparse_binner(unittest.TestCase):

    def test_sum_moduli(self):
        b = binner.SparseBinner_1_1()

        b.accumulate(4j + 2, 2.5, 2)
        b.accumulate(10, 4.0, 0)
        b.accumulate(1j, 10.0, 0)

        self.assertEqual(b.sum_moduli(), sqrt(20.) + 11.)

if __name__ == '__main__':
    unittest.main()
