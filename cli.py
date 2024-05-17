# cli.py

import click
import pandas as pd
import numpy as np
from utils import detectDelimiter
from plot import plot_volcano


@click.command()
@click.option(
    "-f", "--file_path", type=click.Path(exists=True), help="Path to the CSV file"
)
@click.option(
    "-t",
    "--threshold",
    default=7,
    help="Threshold value for significant Fold Change points",
)
@click.option(
    "-o",
    "--output",
    type=str,
    default=".",
    help="Path to save the chart and names file",
)
@click.option(
    "-p",
    "--prefix",
    type=str,
    default="kirekhar",
    help="prefix",
)
def main(file_path, threshold, output, prefix):
    if file_path is None:
        click.echo(
            "Please provide the path to the CSV file using -f or --file_path option."
        )
        return

    delimiter = detectDelimiter(file_path)
    if not delimiter:
        click.echo(
            "Delimiter could not be detected. Please ensure the CSV file has a proper delimiter."
        )
        return

    df = pd.read_csv(file_path, sep=delimiter)

    logFC = df["logFC"]
    P_Value = df["P.Value"]
    neg_log_P_Value = -np.log10(P_Value)
    FC_significant_threshold = threshold

    colors = np.where(
        logFC > FC_significant_threshold,
        "red",
        np.where(logFC < -FC_significant_threshold, "blue", "gray"),
    )

    filtered_lines = df[(abs(logFC) > FC_significant_threshold)]
    with open(f"{output}/{prefix}csv_to_bed.bed", "w") as f:
        for index, row in filtered_lines.iterrows():
            f.write(f"{row['isoform_name']}\n")

    plot_volcano(logFC, neg_log_P_Value, colors, output, prefix)


if __name__ == "__main__":
    main()
