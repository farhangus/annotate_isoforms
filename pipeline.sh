#!/bin/bash

# Function to display usage information
display_annotate_pipeline_usage() {
    # Display options in red color
    echo -e "\e[31mOptions:\e[0m"
    echo -e "\e[31m  -h, --help             Display this help message\e[0m"
    echo -e "\e[31m  -o, --output           Output file path (default: current directory)\e[0m"
    echo -e "\e[31m  -f, --fraction         Similarity fraction (default: 0.9)\e[0m"
    echo -e "\e[31m  -m, --margin           margin base pairs (default: 0)\e[0m"
    echo -e "\e[31m  -f, --fdr              FDR Threshold (default: 0.3) \e[0m"
    echo -e "\e[31m  -c, --FC-threshold     Fold Change (FC) threshold (default: 7)\e[0m"
    echo -e "\e[31m  -p, --prefix           Prefix for output files (default: \"annotated_\")\e[0m"
    echo -e "\e[31m  -i, --input-file       <isoforms_list.names | isoforms_table.csv> Input isoforms list file\e[0m"
    echo -e "\e[31m  -b, --bed              <provided_bedfile> Provided BED file\e[0m"
    echo -e "\e[31m  -r, --reference        <REF> Reference FASTA file\e[0m"
    
    # Display sample command in blue color
    echo -e "\e[34mSample command to run:\e[0m"
    echo -e "\e[34m./pipeline.sh -i list_isoforms.names -b Trinity.bed -r Reference.bed -o annotation_results -p Trinity_\e[0m"
    
    # Additional note
    echo -e "\e[0m** If the provided file is a CSV, it must include columns titled logFC, FDR, and isoform_name.\e[0m"
    echo -e "\e[0m** After adding the margin, the start and end positions of the isoforms will be modified in the output based on the given margin.\e[0m"

}


# Initialize variables
ISOFORM_INPUT_FILE=""
PROVIDED_BEDFILE=""
REFERENCE=""
OUTPUT_FILE="$(pwd)"
PREFIX="annotated_"
FRACTION=0.9
FOLD_CHANGE_THRESHOLD=7
FDR_THRESHOLD=.3
MARGIN=0
CSV_FLAG=0
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
        -c|--FC-threshold)
            FOLD_CHANGE_THRESHOLD=$2
            shift 2
            ;;
        -r|--reference)
            REFERENCE=$2
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE=$2
            rm -rf ${OUTPUT_FILE}
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
        -m|--margin)
            MARGIN=$2
            shift 2
            ;;
        -d|--fdr)
            FDR_THRESHOLD=$2
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

# Check if input file exists
if [[ ! -f "${ISOFORM_INPUT_FILE}" ]]; then
    echo "Error: Input file '${ISOFORM_INPUT_FILE}' does not exist."
    display_annotate_pipeline_usage
    exit 1
fi

# Process the input file if it's a CSV file
if [[ "${ISOFORM_INPUT_FILE##*.}" == "csv" ]]; then
    ROW_ISOFORM_INPUT_FILE=$ISOFORM_INPUT_FILE
    CSV_FLAG=1
    python3 cli.py -f "${ISOFORM_INPUT_FILE}" -t "${FOLD_CHANGE_THRESHOLD}" -o ${OUTPUT_FILE} -p ${PREFIX} -d ${FDR_THRESHOLD}
    ISOFORM_INPUT_FILE="${OUTPUT_FILE}/${PREFIX}csv_to_bed.bed"
fi

# Check if required arguments are provided
if [[ -z $ISOFORM_INPUT_FILE || -z $PROVIDED_BEDFILE || -z $REFERENCE ]]; then
    echo "Missing required arguments."
    display_annotate_pipeline_usage
    exit 1
fi

# Temporary file for isoform locations
tmp_isoforms_locations=$(mktemp)
new_isoforms_locations=$(mktemp)
isoforms_locations="${OUTPUT_FILE}/${PREFIX}isoforms_locations.txt"

# Remove empty lines from input file and write to temp file
grep -v '^\s*$' "$ISOFORM_INPUT_FILE" > ${tmp_isoforms_locations}
rm -f $isoforms_locations

# Iterate through each isoform in the input file
while IFS= read -r isoform; do
    # Check if the isoform exists in the provided BED file
    if grep -q -w "$isoform" "$PROVIDED_BEDFILE"; then
        # Append matching lines from BED file to isoform locations
        grep -w "$isoform" "$PROVIDED_BEDFILE" >> $isoforms_locations
    else 
        # If not found, add a placeholder entry
        echo -e "_\t _\t _\t ${isoform}" >> $isoforms_locations
    fi
done < ${tmp_isoforms_locations}

# Separate found and not found isoforms
grep '^_*_*_' $isoforms_locations > "${OUTPUT_FILE}/${PREFIX}Not_found_isoforms.txt"
ISOFORMS_NOT_FOUND=$(grep -c '^_*_*_' $isoforms_locations)

# Remove not found isoforms from the locations file
grep -v '^_*_*_' ${isoforms_locations} > ${tmp_isoforms_locations}
cat ${tmp_isoforms_locations} > ${isoforms_locations}
sed '/^$/d' ${isoforms_locations} > ${tmp_isoforms_locations}
awk -v margin="$MARGIN" '{start = ($2 - margin < 0) ? 0 : $2 - margin; print $1, start, $3 + margin, $4}' ${tmp_isoforms_locations}  > ${isoforms_locations}
sed -i 's/ \+/	/g' ${isoforms_locations}


# Log and intersect with reference
echo -e "\e[31mTrinity Isoforms found but not found in RefSeq \e[34m Fraction: ${FRACTION}\e[0m" > "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
bedtools intersect -a "$isoforms_locations" -b mm39_genes.bed -f $FRACTION -wa -wb -v >> "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
echo -e "\e[31mIsoforms intersection between Trinity and RefSeq \e[34m Fraction: ${FRACTION}\e[0m" >> "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
bedtools intersect -a "$isoforms_locations" -b mm39_genes.bed -f $FRACTION -wa -wb >> "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt"
awk '!seen[$1,$2,$3,$4]++' "${OUTPUT_FILE}/${PREFIX}log_all_matches.txt" > "${OUTPUT_FILE}/${PREFIX}isoform_result.txt"

# Generate stats
echo -e "\e[31m***************************************************************************************\e[0m" > "${OUTPUT_FILE}/${PREFIX}STATS.txt"
NUMBER_OF_ISOFORMS=$(grep -vc '^\s*$' ${ISOFORM_INPUT_FILE})
awk '!seen[$4]++' "${OUTPUT_FILE}/${PREFIX}isoform_result.txt" > "${OUTPUT_FILE}/${PREFIX}uniq_isoform_result.txt"
echo "Fraction Similarity: ${FRACTION}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Total Number of Isoforms: ${NUMBER_OF_ISOFORMS}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Isoforms Not Found: ${ISOFORMS_NOT_FOUND}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
awk '{print $4}' "${OUTPUT_FILE}/${PREFIX}isoform_result.txt" | sort | uniq -c | awk '$1 > 1' > "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms_names.txt"
awk '{print $2}' "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms_names.txt" | sort | uniq | grep -Ff - "${OUTPUT_FILE}/${PREFIX}isoform_result.txt" > "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms.txt"
rm -f "${OUTPUT_FILE}/${PREFIX}duplicated_isoforms_names.txt"

# Calculate additional stats
ISOFORMS_TRINITY=$(awk '/Trinity Isoforms found/{flag=1; next} /Isoforms intersection between Trinity and RefSeq/{ exit} {if(flag) {count++}} END {print count}' "${OUTPUT_FILE}/${PREFIX}isoform_result.txt")
ISOFORMS_TRINITY_UNIQ=$(awk '/Trinity Isoforms found/{flag=1; next} /Isoforms intersection between Trinity and RefSeq/{ exit} {if(flag) {count++}} END {print count}' "${OUTPUT_FILE}/${PREFIX}uniq_isoform_result.txt")
total_isoform_found=$(wc -l < "${OUTPUT_FILE}/${PREFIX}isoform_result.txt")
total_isoform_found_uniq=$(wc -l < "${OUTPUT_FILE}/${PREFIX}uniq_isoform_result.txt")
echo "Number of Isoforms Exclusively Reported in Trinity: ${ISOFORMS_TRINITY}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Number of Shared Isoforms: $((total_isoform_found - ISOFORMS_TRINITY - 2))" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Number of Unique Isoforms Exclusively Reported in Trinity: ${ISOFORMS_TRINITY_UNIQ}" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"
echo "Number of Shared Unique Isoforms: $((total_isoform_found_uniq - ISOFORMS_TRINITY_UNIQ - 2))" >> "${OUTPUT_FILE}/${PREFIX}STATS.txt"

# genarte the results of exlcusve trinity isoforms with down and up regualted amounts
awk '/Trinity Isoforms found/{flag=1; next} /Isoforms intersection between Trinity and RefSeq/{flag=0} flag' \
    "${OUTPUT_FILE}/${PREFIX}uniq_isoform_result.txt" > "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only.txt"

if [ CSV_FLAG==1 ]
then
    if [ -s "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only.txt" ]; then
        while IFS= read -r line; do
            echo -ne "$line     " >> "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only_up_down_regulated.txt"
            isoform=$(awk '{print$4}' <<< $line) 
            # echo $isoform 
            grep -w $isoform  ${ROW_ISOFORM_INPUT_FILE} >> "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only_up_down_regulated.txt"
        done < "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only.txt"
    else
        echo "The file is empty"
    fi
sed -i 's/,/\t/g' "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only_up_down_regulated.txt"
cat "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only_up_down_regulated.txt"
awk '{print$6}' "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only_up_down_regulated.txt"
while read -r value; do
    # Check if the value is numeric
    if [[ "$value" =~ ^-?[0-9]*\.?[0-9]+$ ]]; then
        # Check if the value is positive
        if (( $(echo "$value > 0" | bc -l) )); then
            (( positive++ ))
        # Check if the value is negative
        elif (( $(echo "$value < 0" | bc -l) )); then
            (( negative++ ))
        fi
    fi
done <<< "$(awk '{print$6}' "${OUTPUT_FILE}/${PREFIX}uniq_isoforms_trinty_only_up_down_regulated.txt")"
echo "upregulated: $positive, downregulated: $negative"

fi