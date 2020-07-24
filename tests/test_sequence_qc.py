#!/usr/bin/env python

import os
import pytest
from pytest import approx

from sequence_qc.noise import calculate_noise, OUTPUT_NOISE_FILENAME, OUTPUT_PILEUP_NAME


def test_calculate_noise():
    """
    Test noise calculation from pysamstats

    :return:
    """
    noise = calculate_noise(
        'test_data/ref_nochr.fa',
        'test_data/SeraCare_0-5.bam',
        'test_data/test.bed',
        0.2,
        output_prefix='test_'
    )
    assert noise == approx(0.0012269938650291694, rel=1e-6)

    for filename in [
            'test_' + OUTPUT_PILEUP_NAME,
            'test_' + OUTPUT_NOISE_FILENAME,
            'test_noise_acgt.tsv',
            'test_noise_del.tsv',
            'test_noise_n.tsv',
    ]:
        assert os.path.exists(filename)
        os.unlink(filename)


if __name__ == '__main__':
    pytest.main()
