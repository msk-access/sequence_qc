#!/usr/bin/env python

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
        noise = sequence_qc.calculate_noise(
            'test_data/ref.fa',
            'test_data/SeraCare_0-5.bam',
            'test_data/test.bed',
            0.2
        )
        assert noise == 0.0012239902079405065
