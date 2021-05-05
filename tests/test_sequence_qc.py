#!/usr/bin/env python

import os
import pytest
from pytest import approx
import pandas as pd

from sequence_qc.noise import calculate_noise, OUTPUT_NOISE_FILENAME, OUTPUT_PILEUP_NAME
from sequence_qc import plots
from sequence_qc.noise_by_tlen import get_fragment_size_for_sample

CUR_DIR = os.path.dirname(os.path.abspath(__file__))


def test_calculate_noise():
    """
    Test noise calculation from pysamstats

    :return:
    """
    noise = calculate_noise(
        os.path.join(CUR_DIR, 'test_data/ref_nochr.fa'),
        os.path.join(CUR_DIR, 'test_data/SeraCare_0-5_tmp.bam'),
        os.path.join(CUR_DIR, 'test_data/test.bed'),
        0.2,
        sample_id='test_'
    )
    # noise = calculate_noise(
    #     os.path.join(CUR_DIR, '/Users/johnsoni/Desktop/test_bam/Homo_sapiens_assembly19.fasta'),
    #     os.path.join(CUR_DIR, '/Users/johnsoni/Desktop/test_bam/C-XJ1562-L012-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'),
    #     os.path.join(CUR_DIR, '/Users/johnsoni/Desktop/test_bam/MSK-ACCESS-v1_0-probe-A_no_msi_sorted_deduped.bed'),
    #     0.02,
    #     sample_id='C-XJ1562-L012-d_'
    # )
    assert noise == approx(0.0048899755501162715, rel=1e-6)

    for filename in [
            'test_' + OUTPUT_PILEUP_NAME,
            'test_' + OUTPUT_NOISE_FILENAME,
            'test_noise_acgt.tsv',
            'test_noise_del.tsv',
            'test_noise_n.tsv',
    ]:
        assert os.path.exists(filename)
        os.unlink(filename)


def test_noise_by_tlen():
    """
    """
    noisy_positions = pd.read_csv(os.path.join(CUR_DIR, 'test_data/SeraCare_noise_positions.tsv'), sep='\t')
    get_fragment_size_for_sample(
        'test',
        os.path.join(CUR_DIR, 'test_data/SeraCare_0-5_tmp.bam'),
        'test',
        noisy_positions, 0, 500
    )


def test_noisy_positions_plot():
    """
    Test HTML plot from plotly is produced

    :return:
    """
    noise_df = pd.read_csv(os.path.join(CUR_DIR, 'test_data/test_noise_positions.tsv'), sep='\t')
    plots.plot_noisy_positions(noise_df)


def test_n_counts_plot():
    """
    Test HTML plot for N counts

    :return:
    """
    noise_df = pd.read_csv(os.path.join(CUR_DIR, 'test_data/test_noise_positions.tsv'), sep='\t')
    plots.plot_n_counts(noise_df)


def test_all_plots():
    """
    Test combined HTML plot

    :return:
    """
    noise_df = pd.read_csv(os.path.join(CUR_DIR, 'test_data/test_noise_positions.tsv'), sep='\t')
    noise_by_substitution = pd.read_csv(
        os.path.join(CUR_DIR, 'test_data/test_noise_by_substitution.tsv'), sep='\t')
    filename = os.path.join(CUR_DIR, 'test_data/test_noise_by_tlen.tsv')

    noisy_tlen_df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["Sample", "Type", "read_id", "Var", "Size", "Chr", "Pos", "geno_not_geno"],
        dtype="object",
    )

    plots.all_plots(noise_df, noise_df, noise_by_substitution, noisy_tlen_df)
    assert os.path.exists('_noise.html')


if __name__ == '__main__':
    pytest.main()
