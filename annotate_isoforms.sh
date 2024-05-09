#!/bin/bash

while IFS= read -r line; do     isoform=$(echo "$line");grep  "$isoform" sorted_Trinity_all_samples_gmap.bed ;  done < $1
