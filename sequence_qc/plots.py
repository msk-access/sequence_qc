import pandas as pd
import plotly
import plotly.express as px

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sequence_qc.plot_noise_by_tlen import create_noisy_tlen_plot

# from sequence_qc.noise import NOISE_FRACTION


def all_plots(pileup_df: pd.DataFrame, noisy_positions: pd.DataFrame, st_df: pd.DataFrame, noisy_tlen_df: pd.DataFrame,
              sample_id: str = '') -> None:
    """
    Create all plots in a single HTML report

    :param pileup_df: pd.DataFrame - All positions from bed file as data frame
    :param noisy_positions: pd.DataFrame - Noisy positions data frame
    :param st_df: pd.DataFrame - Substitution types data frame
    :param noisy_tlen_df: pd.DataFrame -
    :param sample_id:
    :return:
    """
    with open(sample_id + '_noise.html', 'w') as f:
        f.write('<h1 style=\'font-family: sans-serif\'>Noise Report for sample {}</h1>'.format(sample_id))

        fig = plot_noise_by_substitution(st_df)
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        fig = plot_noisy_positions(noisy_positions)
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        fig = create_noisy_tlen_plot(noisy_tlen_df)
        if fig:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        fig = plot_n_counts(pileup_df)
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))


def plot_noise_by_substitution(st_df: pd.DataFrame) -> plotly.graph_objects.Figure:
    """
    Barplot for noise fraction of each substitution type

    :param st_df: pd.DataFrame - Data frame with rows for each substitution type, and column 'NOISE_FRACTION'
    :return:
    """
    title = 'Noise By Substitution'
    fig = px.bar(
        x=st_df.index,
        # Making units easier to display by raising to 10^6
        y=st_df['noise_fraction'] * (10**6),
        title=title,
        labels={'x': 'Substitution', 'y': 'Alt Count / (Alt Count + Ref Count) x 10^6'}
    )
    fig.update_layout(
        yaxis=dict(
            showexponent='none',
            exponentformat='none'
        )
    )
    return fig


def plot_noisy_positions(noisy_pileup_df: pd.DataFrame) -> plotly.graph_objects.Figure:
    """
    Barplot and violin plot of positions with most noise, as defined by calculate_noise module

    :param noisy_pileup_df:
    :return:
    """
    noisy_pileup_df = noisy_pileup_df.sort_values('noise_acgt', ascending=False)
    noisy_pileup_df['chrom_pos'] = noisy_pileup_df['chrom'] + ':' + noisy_pileup_df['pos'].astype(str)
    bar_title = 'Top 100 Noisy Positions'
    box_title = 'All positions'

    fig = make_subplots(
        rows=1,
        cols=4,
        specs=[[{"colspan": 3}, None, None, {}]],
        subplot_titles=(bar_title, box_title),
    )
    fig.update_layout(showlegend=False)
    fig.update_yaxes(title_text="Alt Count / Total Count", range=[0, 0.2], row=1, col=1)
    fig.update_xaxes(title_text="Genomic Position", row=1, col=1, showticklabels=False)

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
    fig.update_xaxes(title_text="", row=1, col=4, showticklabels=False)
    return fig


def plot_noise_by_tlen(avg_tlen_noise: pd.DataFrame) -> px.bar:
    """
    Barplot of average template length at noisy sites
    """
    mean_tlen_fwd = avg_tlen_noise['mean_tlen_fwd'].value_counts()
    mean_tlen_rev = avg_tlen_noise['mean_tlen_rev'].value_counts()
    mean_tlen = pd.concat([mean_tlen_fwd, mean_tlen_rev])
    title = 'Mean template length for positions with noise (calculated separately for fwd/rev reads)'
    fig = px.bar(
        x=mean_tlen.index,
        y=mean_tlen,
        title=title,
        labels={'x': 'TLEN', 'y': 'Number of noisy positions'}
    )
    fig.update_xaxes(range=[-200, 200])
    return fig


def plot_n_counts(pileup_df: pd.DataFrame) -> px.bar:
    """
    Barplot of number of sites with each discrete N count

    :param pileup_df:
    :return:
    """
    n_counts = pileup_df['N'].value_counts()
    title = 'Positions with each N count'
    fig = px.bar(
        x=n_counts.index,
        y=n_counts,
        title=title,
        labels={'x': 'N count', 'y': 'Number of positions'}
    )
    fig.update_xaxes(range=[0, 50])
    return fig
