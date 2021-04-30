import numpy as np
import pandas as pd

import plotly.figure_factory as ff


SUBSTITUTION_TYPES_COMBINED = [
    [["G>T", "C>A"], "C>A"],
    [["C>G", "G>C"], "C>G"],
    [["G>A", "C>T"], "C>T"],
    [["T>A", "A>T"], "T>A"],
    [["A>G", "T>C"], "T>C"],
    [["T>G", "A>C"], "T>G"],
]

COLORS = ["#556278", "#C1292E", "#F2A535"]


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


def plot_data(frag_size_select_df):
    """
    Select data for creating plot

    :param: frag_size_select_df - pd.DataFrame
    :param:
    """
    frag_size_select_df["Size"] = pd.to_numeric(
        frag_size_select_df["Size"], downcast="integer"
    )
    frag = frag_size_select_df[(frag_size_select_df["Size"] > 0)]
    # write_summary(frag, outfile)
    fig = plot_data_type(frag)
    return fig


def plot_data_type(frag):
    """
    Create the plot

    :param: frag pd.DataFrame -
    :param:
    """
    data = []
    labels = []
    for st_pair, st in SUBSTITUTION_TYPES_COMBINED:
        st_sizes = frag[frag['Var'].isin(st_pair)]['Size']
        if len(st_sizes) > 1:
            data.append(st_sizes)
            labels.append(st)
    geno_series = frag[frag['Var'] == 'GENOTYPE']['Size']
    if len(geno_series) > 1:
        data.append(geno_series)
        labels.append('Genotype')
    n_series = frag[frag['Var'] == 'N']['Size']
    if len(n_series) > 1:
        data.append(n_series)
        labels.append('N')

    try:
        fig = ff.create_distplot(
            data,
            labels,
            # colors=COLORS,
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
    frag_size_select_df = noisy_tlen_df[['Var', 'Size', 'Chr', 'Pos']]
    fig = plot_data(frag_size_select_df)
    return fig
