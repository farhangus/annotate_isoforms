# Annotate_Trinity_isoforms

Annotate_Trinity_isoforms is a mixture of Bash and Python packages that, through using Python libraries, along with the `bedtools`, allows us to annotate all detected isoforms by the Trinity pipeline. It prepares the list of exclusive isoforms detected by Trinity not reported in the reference genome and also prepares the list of shared isoforms between Trinity and the reference genome.

## Prerequisites

- Python 3.x
- `click` library
- `pandas` library
- `numpy` library
- `bedtools` installed and added to your system's PATH

## Sample Usage

```sh
./pipeline.sh -i DEG_FULL.csv -o tmp/sample_output_folder --bed sorted_Trinity_all_samples_gmap.bed --reference mm39_genes.bed --fraction 0.9 --FC-threshold 6 --fdr 0.01

