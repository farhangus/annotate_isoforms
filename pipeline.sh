#!/bin/bash

# Function to display usage information
display_annotate_pipeline_usage() {
# Red color
echo -e "\e[31mOptions:\e[0m"
echo -e "\e[31m  -h,--help             Display this help message\e[0m"
echo -e "\e[31m  -o,--output           output_file_path                            default=current_directory\e[0m"
echo -e "\e[31m  -f,--fraction         Similarity fraction                         default=0.9\e[0m"
echo -e "\e[31m  -c,--FC-threshhold    FC threshhold                               default=7\e[0m"
echo -e "\e[31m  -p,--prefix           prefix                                      default=\"annotated_\"\e[0m"
echo -e "\e[31m  -i,--input-file       <isoforms_list.names | isoforms_table.csv>  Input isoforms_list file\e[0m"
echo -e "\e[31m  -b,--bed              <provided_bedfile>                          Provided BED file\e[0m"
echo -e "\e[31m  -r,--reference        <REF>                                       Reference FASTA file\e[0m"
# Blue color
echo -e "\e[34mSample command for run:\e[0m"
echo -e "\e[34m./pipeline.sh -i list_isoforms.names -b Trinity.bed -r Reference.bed -o annotation_results -p Trinity_\e[0m"
# Reset color
echo -e "\e[0m"
}

# Initialize variables
ISOFORM_INPUT_FILE=""
PROVIDED_BEDFILE=""
REFERENCE=""
OUTPUT_FILE="$(pwd)"
PREFIX="annotated_"
TOTAL_NUMBER_OF_ISOFORMS=0
FOUND_ISOFORMS=0
FRACTION=.9
NUMBER_OF_ISOFORMS=""
ISOFORMS_NOT_FOUND=""
ISOFORMS_TRINITY=""
FOLD_CHANGE_THRESHHOLD="7"

# Parse command-line options
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-file)
            ISOFORM_INPUT_FILE=$2
            shift 2
            ;;
        -b|--bed)
            PROVIDED_BEDFILE=$2
            shift 2
            ;;
         -c|--FC-threshhold )
            FOLD_CHANGE_THRESHHOLD=$2
            shift 2
            ;;
        -r|--reference)
            REFERENCE=$2
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE=$2
            rm -r ${OUTPUT_FILE}
            mkdir -p ${OUTPUT_FILE}
            shift 2
            ;;
        -p|--prefix)
            PREFIX=$2
            shift 2
            ;;
        -f|--fraction)
            FRACTION=$2
            shift 2
            ;;
        -h|--help)
            display_annotate_pipeline_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            display_annotate_pipeline_usage
            exit 1
            ;;
    esac
done
if [[ ! -f "${ISOFORM_INPUT_FILE}" ]]; then
    echo "Error: Input file "${ISOFORM_INPUT_FILE}" does not exist."
    display_annotate_pipeline_usage
    exit 1
fi
# Check if the file is a CSV file based on its extension

if [[ "${ISOFORM_INPUT_FILE##*.}" == "csv" ]]; then
    python3 report_isoforms_names.py -f "${ISOFORM_INPUT_FILE}" -t "${FOLD_CHANGE_THRESHHOLD}"  -o ${OUTPUT_FILE}
    ISOFORM_INPUT_FILE="${OUTPUT_FILE}/csv_to_bed.bed"
fi
# Check if required arguments are provided
if [[ -z $ISOFORM_INPUT_FILE || -z $PROVIDED_BEDFILE || -z $REFERENCE ]]
then
    echo "Missing required arguments."
    display_annotate_pipeline_usage
    exit 1
fi
tmp_isoforms_locations=$(mktemp)
isoforms_locations="${OUTPUT_FILE}/${PREFIX}isoforms_locations.txt"
grep -v '^\s*$' "$ISOFORM_INPUT_FILE" >${tmp_isoforms_locations}
rm -f $isoforms_locations
# Iterate through each line in the ISOFORM_INPUT_FILE file
while IFS= read -r isoform; do
    # Check if the isoform exists in the provided BED file
    if grep -q "$isoform" "$PROVIDED_BEDFILE"; then
        # Get the lines from the BED file matching the isoform
        grep "$isoform" "$PROVIDED_BEDFILE" >> $isoforms_locations
    else 
        echo -e "_\t _\t _\t ${isoform}" >> $isoforms_locations
    fi
done < ${tmp_isoforms_locations}
grep '^_*_*_' $isoforms_locations > "${OUTPUT_FILE}/${PREFIX}Not_found_isoforms.txt"
ISOFORMS_NOT_FOUND=$(grep -c '^_*_*_' $isoforms_locations)
 
grep -v '^_*_*_' ${isoforms_locations} > ${tmp_isoforms_locations}
cat ${tmp_isoforms_locations} > ${isoforms_locations}

echo -e "\e[31mTrinity Isoforms found  but not found in RefSeq \e[34m Fraction: ${FRACTION}\e[0m" > "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
bedtools intersect -a "$isoforms_locations" -b mm39_genes.bed -f $FRACTION  -wa -wb  -v >> "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
echo -e "\e[31mIsoforms intersection between Trinity and RefSeq \e[34m Fraction: ${FRACTION}\e[0m" >> "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
bedtools intersect -a "$isoforms_locations" -b mm39_genes.bed -f $FRACTION  -wa -wb  >> "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
awk '!seen[$1,$2,$3]++' "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"  > "${OUTPUT_FILE}/${PREFIX}isoform_result.txt"
echo -e "\e[31m***************************************************************************************\e[0m" > "${OUTPUT_FILE}/${PREFIX}STATS.txt"
NUMBER_OF_ISOFORMS="$(grep -v '^\s*$' ${ISOFORM_INPUT_FILE} | wc -l  | grep -o '[0-9]*' )"
awk '!seen[$4]++' "${OUTPUT_FILE}/${PREFIX}isoform_result.txt" > "${OUTPUT_FILE}/${PREFIX}uniq_isoform_result.txt"
echo "Fraction Similarity: ${FRACTION}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Total Number of Isoforms: ${NUMBER_OF_ISOFORMS}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Isoforms Not Found: ${ISOFORMS_NOT_FOUND}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
awk '{print $4}' "${OUTPUT_FILE}/${PREFIX}isoform_result.txt" | sort | uniq -c | awk '$1 > 1' > "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms_names.txt"
awk '{print$2}' "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms_names.txt" | sort | uniq | grep -Ff - "${OUTPUT_FILE}/${PREFIX}isoform_result.txt" > "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms.txt"
rm "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms_names.txt"

ISOFORMS_TRINITY=$(awk '/Trinity Isoforms found /{flag=1; next} /Isoforms intersection between Trinity and RefSeq/{ exit} {if(flag) {count++}} END {print count}' "${OUTPUT_FILE}/${PREFIX}isoform_result.txt")
ISOFORMS_TRINITY_UNIQ=$(awk '/Trinity Isoforms found /{flag=1; next} /Isoforms intersection between Trinity and RefSeq/{ exit} {if(flag) {count++}} END {print count}' "${OUTPUT_FILE}/${PREFIX}uniq_isoform_result.txt")
echo "Number of Isoforms Exclusively Reported in Trinity: ${ISOFORMS_TRINITY}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
total_isoform_found=$(wc -l "${OUTPUT_FILE}/${PREFIX}isoform_result.txt"| grep -o '[0-9]*')
total_isoform_found_uniq=$(wc -l "${OUTPUT_FILE}/${PREFIX}uniq_isoform_result.txt"| grep -o '[0-9]*')
echo "Number of Shared Isoforms: $((total_isoform_found -ISOFORMS_TRINITY -2 ))" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Number of Unique Isoforms Exclusively Reported in Trinity: ${ISOFORMS_TRINITY_UNIQ}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Number of Shared Unique Isoforms: $((total_isoform_found_uniq - ISOFORMS_TRINITY_UNIQ -2 ))" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"

