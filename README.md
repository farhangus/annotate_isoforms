# Annotate_Trinity_isoforms

Annotate_Trinity_isoforms is a mixture of Bash and Python packages that, through the use of Python libraries along with the `bedtools`, allows us to annotate all detected isoforms by the Trinity pipeline. It prepares the list of exclusive isoforms detected by Trinity not reported in the reference genome and also prepares the list of shared isoforms between Trinity and the reference genome. The pipeline generates the following files:

- **annotated_uniq_isoforms_trinty_only_up_down_regulated.txt**
  - Description: Contains uniquely annotated isoforms detected by Trinity, with up/down-regulated information. It also reports the total number of down-regulated versus upregulated isoforms.

- **annotated_uniq_isoforms_trinty_only.txt**
  - Description: Contains uniquely annotated isoforms detected by Trinity only.

- **annotated_STATS.txt**
  - Description: Contains statistics and summary information from the annotation process. It provides the following information:
    - Fraction Similarity
    - Total Number of Isoforms
    - Isoforms Not Found
    - Number of Isoforms Exclusively Reported in Trinity
    - Number of Shared Isoforms
    - Number of Unique Isoforms Exclusively Reported in Trinity
    - Number of Shared Unique Isoforms

- **annotated_duplicated_isoforms.txt**
  - Description: Contains duplicated annotated isoforms, which are reported in two different locations.

- **annotated_uniq_isoform_result.txt**
  - Description: Contains the result of annotated unique isoforms.

- **annotated_isoform_result.txt**
  - Description: Contains the result of annotated isoforms.

- **annotated_log_all_matches.txt**
  - Description: Contains logs of all matches during the annotation process.

- **annotated_isoforms_locations.txt**
  - Description: Contains the locations of annotated isoforms.

- **annotated_Not_found_isoforms.txt**
  - Description: Contains the isoforms that were not found during the annotation process.

- **annotated_volcano_grouped.png**
  - Description: PNG image showing a grouped volcano plot.

## Prerequisites

- Python 3.x
- `click` library
- `pandas` library
- `numpy` library
- `bedtools` installed and added to your system's PATH

## Sample Usage

```sh
For the provided CSV file, the columns named logFC, FDR, and isoform_name are required. If your CSV file has different names for these columns, please manipulate them accordingly.
./pipeline.sh --input-file DEG_FULL.csv --output tmp_sample_output/ --bed sorted_Trinity_all_samples_gmap.bed --reference mm39_genes.bed --similarity-fraction 0.90 --FC-threshold 6 --fdr 0.01

