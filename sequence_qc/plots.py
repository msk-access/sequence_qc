import pandas as pd
import plotly.express as px


TOP_NOISE_PLOT = 'noisy_positions.html'


def plot_top_noisy_positions(noisy_pileup_df: pd.DataFrame, output_prefix: str = '') -> None:
    """
    Barplot of positions with most noise, as defined by calculate_noise module

    :param pileup_df:
    :return:
    """
    noisy_pileup_df = noisy_pileup_df.sort_values('noise_acgt')
    noisy_pileup_df['chrom_pos'] = noisy_pileup_df['chrom'] + ':' + noisy_pileup_df['pos'].astype(str)
    fig = px.bar(noisy_pileup_df, x='chrom_pos', y='noise_acgt')
    fig.write_html(output_prefix + TOP_NOISE_PLOT)
