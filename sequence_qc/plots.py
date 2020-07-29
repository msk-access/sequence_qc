import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


TOP_NOISE_PLOT = 'noisy_positions.pdf'


def plot_top_noisy_positions(noisy_pileup_df: pd.DataFrame) -> None:
    """
    Barplot of positions with most noise, as defined by calculate_noise module

    :param pileup_df:
    :return:
    """
    plt.clf()
    plt.figure(figsize=(10, 5))
    pos = np.arange(len(noisy_pileup_df))
    noise = noisy_pileup_df['noise_acgt']
    x_lab = noisy_pileup_df['chrom'] + ':' + noisy_pileup_df['pos'].astype(str)
    plt.bar(pos, noise, align="center", color="black")
    plt.xticks(ticks=pos, labels=x_lab, rotation=90, ha="center")
    plt.title('Top Noisy Positions')
    plt.savefig(TOP_NOISE_PLOT)
