#!/usr/bin/env python

"""Tests for `sequence_qc` package."""


import unittest

from sequence_qc import sequence_qc


class TestSequence_qc(unittest.TestCase):
    """Tests for `sequence_qc` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_calculate_noise(self):
        """
        Test noise calculation

        :return:
        """
        noise = sequence_qc.calculate_noise('test_data/SeraCare_0-5.bam', 'test_data/test.bed', 0.002)
        assert noise == 0.00001
