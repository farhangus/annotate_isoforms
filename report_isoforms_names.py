#!/usr/bin/python3

import click
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define a command-line interface using Click
@click.command()
@click.option('-f', '--file_path', type=click.Path(exists=True), help='Path to the CSV file')
@click.option('-t', '--threshold', default=7, help='Threshold value for significant Fold Change points')
@click.option('-o', '--output', type=str, default='.', help='Path to save the chart and names file')
def main(file_path, threshold, output):
    # Check if the file path is provided
    if file_path is None:
        click.echo("Please provide the path to the CSV file using -f or --file_path option.")
        return
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path,sep='\s+')
    #df = pd.read_csv(file_path, sep=",", engine='python')
    # Extract columns from DataFrame
    logFC = df['logFC']
    P_Value = df['P.Value']
    print(logFC)
    # Calculate -log10(P.Value)
    neg_log_P_Value = -np.log10(P_Value)

    # Set threshold for significant points
    FC_significant_threshold = threshold

    # Define colors based on conditions
    colors = np.where(logFC > FC_significant_threshold, 'red', np.where(logFC < -FC_significant_threshold, 'blue', 'gray'))

    # Define legend labels
    legend_labels = {
        'gray': 'Non-Significant',
        'red': 'Positive FC',
        'blue': 'Negative FC'
    }

    # Write significant points to a BED file
    filtered_lines = df[(abs(logFC) > FC_significant_threshold)]
    with open(f'{output}/csv_to_bed.bed', 'w') as f:
        for index, row in filtered_lines.iterrows():
            f.write(f"{row['isoform_name']}\n")

    # Plot volcano plot
    plt.figure(figsize=(8, 6))
    for color in np.unique(colors):
        plt.scatter(logFC[colors == color], neg_log_P_Value[colors == color], color=color, alpha=0.7, label=legend_labels[color])

    # Set plot title and axis labels
    plt.title('Volcano Plot for all samples')
    plt.xlabel('logFC (Fold Change)')
    plt.ylabel('-log10(P.Value)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'{output}/volcano_grouped.png')  # Save the plot as an image file

if __name__ == "__main__":
    main()  # Execute the main function when the script is run
