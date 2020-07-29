import pandas as pd
import plotly.express as px


TOP_NOISE_PLOT = 'noisy_positions.html'


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
