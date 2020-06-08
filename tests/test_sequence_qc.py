#!/usr/bin/env python

import pytest
from pytest import approx

from sequence_qc.noise import calculate_noise


def test_calculate_noise():
    """
    Test noise calculation from pysamstats

    :return:
    """
    noise = calculate_noise(
        'test_data/ref_nochr.fa',
        'test_data/SeraCare_0-5.bam',
        'test_data/test.bed',
        0.2
    )
    assert noise == approx(0.001228501228, rel=1e-6)


if __name__ == '__main__':
    pytest.main()
