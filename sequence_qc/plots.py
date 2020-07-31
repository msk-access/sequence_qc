import pandas as pd
import plotly.express as px

import plotly.graph_objects as go
from plotly.subplots import make_subplots


TOP_NOISE_PLOT = 'noisy_positions.html'
N_COUNTS_PLOT = 'n_counts.html'


def plot_noisy_positions(noisy_pileup_df: pd.DataFrame, sample_id: str = '') -> None:
    """
    Barplot and violin plot of positions with most noise, as defined by calculate_noise module

    :param pileup_df:
    :return:
    """
    noisy_pileup_df = noisy_pileup_df.sort_values('noise_acgt', ascending=False)
    noisy_pileup_df['chrom_pos'] = noisy_pileup_df['chrom'] + ':' + noisy_pileup_df['pos'].astype(str)
    bar_title = 'Top 100 Noisy Positions for sample {}'.format(sample_id)
    box_title = 'All positions'

    fig = make_subplots(
        rows=1,
        cols=4,
        specs=[[{"colspan": 3}, None, None, {}]],
        subplot_titles=(bar_title, box_title),
    )
    fig.update_layout(showlegend=False)
    fig.update_yaxes(title_text="minor allele count / total", range=[0, 0.2], row=1, col=1)
    fig.update_xaxes(title_text="contig:position", row=1, col=1)

    noise_subset = noisy_pileup_df[:100]

    fig.add_trace(
        go.Bar(
            x=noise_subset['chrom_pos'],
            y=noise_subset['noise_acgt'],
            text=noise_subset['minor_allele_count']
        ),
        row=1, col=1
    )

    fig.add_trace(
        go.Violin(
            y=noisy_pileup_df['noise_acgt'],
        ),
        row=1, col=4
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
