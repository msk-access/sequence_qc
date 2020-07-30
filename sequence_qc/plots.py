import pandas as pd
import plotly.express as px


TOP_NOISE_PLOT = 'noisy_positions.html'
N_COUNTS_PLOT = 'n_counts.html'


def plot_top_noisy_positions(noisy_pileup_df: pd.DataFrame, sample_id: str = '') -> None:
    """
    Barplot of positions with most noise, as defined by calculate_noise module

    :param pileup_df:
    :return:
    """
    noisy_pileup_df = noisy_pileup_df.sort_values('noise_acgt', ascending=False)
    noisy_pileup_df['chrom_pos'] = noisy_pileup_df['chrom'] + ':' + noisy_pileup_df['pos'].astype(str)
    title = 'Top 100 Noisy Positions for sample {}'.format(sample_id)
    fig = px.bar(
        noisy_pileup_df[:100],
        x='chrom_pos',
        y='noise_acgt',
        title=title,
        hover_data=['ref', 'A', 'C', 'G', 'T', 'minor_allele_count', 'major_allele_count']
    )
    fig.write_html(sample_id + TOP_NOISE_PLOT)


def plot_n_counts(pileup_df: pd.DataFrame, sample_id: str = '') -> None:
    """
    Barplot of number of sites with each discrete N count

    :param pileup_df:
    :param sample_id:
    :return:
    """

    n_counts = pileup_df['N'].value_counts()
    title = 'Positions with each N count, for sample {}'.format(sample_id)
    fig = px.bar(
        x = n_counts.index,
        y = n_counts,
        title = title,
        labels = {'x': 'N count', 'y': 'Number of positions'}
    )
    fig.write_html(sample_id + N_COUNTS_PLOT)
