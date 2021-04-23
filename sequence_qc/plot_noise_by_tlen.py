import numpy as np
import pandas as pd

import plotly.figure_factory as ff


def write_summary(df, outfile):
    summary = (
        df.groupby(["Sample", "Var"], as_index=False)[["Type"]]
        .count()
        .rename(columns={"Type": "Count"})
    )
    all_summary = (
        df.groupby(["Var"], as_index=False)[["Type"]]
        .count()
        .rename(columns={"Type": "Count"})
    )
    all_summary["Sample"] = "Total"
    summary = pd.concat([all_summary, summary], ignore_index=True)[
        ["Sample", "Var", "Count"]
    ]
    summary.to_csv(outfile.replace(".pdf", ".txt"), sep="\t", index=False)


def plot_data(frag_size_select_df, outfile):
    """
    Select data for creating plot

    :param: frag_size_select_df - pd.DataFrame
    :param:
    """
    frag_size_select_df["Size"] = pd.to_numeric(
        frag_size_select_df["Size"], downcast="integer"
    )
    frag = frag_size_select_df[(frag_size_select_df["Size"] > 0)]
    write_summary(frag, outfile)
    fig = plot_data_type(frag)
    return fig


def plot_data_type(frag):
    """
    Create the plot

    :param: frag pd.DataFrame -
    :param:
    """
    geno_series = frag[frag['Var'] == 'GENOTYPE']['Size']
    noise_series = frag[frag['Var'] == 'NOISE']['Size']
    data = [geno_series, noise_series]
    try:
        fig = ff.create_distplot(
            data,
            ['Genotype', 'Noise'],
            colors=["#556278", "#C1292E"],
            bin_size=.2,
            show_rug=False,
        )
        fig.update_layout(
            title='Fragment size distribution for noisy positions',
            xaxis_title="TLEN",
        )
    except np.linalg.LinAlgError:
        # Not enough data for plot
        return None
    return fig


def create_noisy_tlen_plot(noisy_tlen_df):
    """
    Interface to this module

    :return:
    """
    frag_size_select_df = noisy_tlen_df[["Sample", "Type", "Var", "Size"]]

    fig = plot_data(
        frag_size_select_df,
        'noise_by_tlen.pdf',
    )
    return fig
