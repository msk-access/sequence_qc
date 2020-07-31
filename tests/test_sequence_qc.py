#!/usr/bin/env python

import os
import pytest
from pytest import approx
import pandas as pd

from sequence_qc.noise import calculate_noise, OUTPUT_NOISE_FILENAME, OUTPUT_PILEUP_NAME
from sequence_qc import plots


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
        sample_id='test_'
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


def test_noisy_positions_plot():
    """
    Test HTML plot from plotly is produced

    :return:
    """
    noise_df = pd.read_csv('test_data/test_noise_positions.tsv', sep='\t')
    plots.plot_noisy_positions(noise_df)


def test_n_counts_plot():
    """
    Test HTML plot for N counts

    :return:
    """
    noise_df = pd.read_csv('test_data/test_noise_positions.tsv', sep='\t')
    plots.plot_n_counts(noise_df)


def test_all_plots():
    """
    Test combined HTML plot

    :return:
    """
    noise_df = pd.read_csv('test_data/test_noise_positions.tsv', sep='\t')
    plots.all_plots(noise_df, noise_df)
    assert os.path.exists('noise.html')


if __name__ == '__main__':
    pytest.main()
