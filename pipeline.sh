#!/bin/bash

# Function to display usage information
display_aanotate_pipeline_usage() {
    echo "Options:"
    echo "  -h, --help                          Display this help message"
    echo "  -f, --file      <isoforms_list>     Input isoforms_list file"
    echo "  -b, --bed       <provided_bedfile>  provided_bedfile "
    echo "  -r, --reference <REF>               reference fasta file"
    echo "  -m, --margin    base_pairs_margin   default=10"
    echo "  -o, --output    output_fle_path   output file path"
}

# Initialize variables
ISOFORMS_LIST=""
PROVIDED_BEDFILE=""
REFERENCE=""
OUTPUT_FILE=""
MARGIN=0

# Parse command-line options
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--file)
            ISOFORMS_LIST=$2
            shift 2
            ;;
        -b|--bed)
            PROVIDED_BEDFILE=$2
            shift 2
            ;;
        -r|--reference)
            REFERENCE=$2
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE=$2
            shift 2
            ;;
        -m|--margin)
            MARGIN=$2
            shift 2
            ;;
        -h|--help)
            display_aanotate_pipeline_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            display_aanotate_pipeline_usage
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z $ISOFORMS_LIST || -z $PROVIDED_BEDFILE || -z $REFERENCE ]]
then
    echo "Missing required arguments."
    display_aanotate_pipeline_usage
    exit 1
fi
isoforms_locations=$(mktemp)
rm -f $isoforms_locations
# Iterate through each line in the ISOFORMS_LIST file
while IFS= read -r isoform; do
    # Check if the isoform exists in the provided BED file
    if grep -q "$isoform" "$PROVIDED_BEDFILE"; then
        # Get the lines from the BED file matching the isoform
        grep "$isoform" "$PROVIDED_BEDFILE" >> $isoforms_locations
    else 
        echo -e "\t\t\t$isoform" >> $isoforms_locations
    fi
done < "$ISOFORMS_LIST"
echo $isoforms_locations
while IFS= read -r line; do
    bedfile_chromosome=$(awk '{print $1}' <<< "$line")
    gene_start_position=$(awk '{print $2}' <<< "$line")   
    gene_end_position=$(awk '{print $3}' <<< "$line")   
    trinity_isoform_name=$(awk '{print $4}' <<< "$line")   
    trinity_gene_name=$(sed 's/_[^_]*$//'<<<$trinity_isoform_name)

    if [[ ! -z "$bedfile_chromosome" && ! -z "$gene_start_position" ]]
    then
    start=$((gene_start_position - MARGIN))
    end=$((gene_start_position + MARGIN))
    found_match=false
    # Iterate over the range
    for ((i = start; i <= end; i++)); do
        if grep -q  -E -m 1 "^$bedfile_chromosome\s+$i\s+" "$REFERENCE"
        then
            line=$(grep -E -m 1 "^$bedfile_chromosome\s+$i\s+" "$REFERENCE")
            echo "$bedfile_chromosome" "$gene_start_position", "$gene_end_position" , "$trinity_isoform_name" , "$trinity_gene_name", "$line"
            found_match=true
            break
 
        else
            line=""
        fi
    done
    if ! $found_match
    then
        echo "$bedfile_chromosome" "$gene_start_position", "$gene_end_position" , "$trinity_isoform_name" , "$trinity_gene_name", "$line"
    fi
    fi  
        
    


done < "$isoforms_locations"