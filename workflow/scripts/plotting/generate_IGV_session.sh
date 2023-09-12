#!/usr/bin/env bash

# if [ "$#" -ne 1 ]; then
#     echo "Usage: $0 <path_to_track_directory>"
#     exit 1
# fi

# track_dir=${snakemake_input[splitted_files_dir]}
track_dir="$1"
# output_file=${snakemake_output[igv_session]}
output_file="$2"

# Check if the directory exists
if [ ! -d "$track_dir" ]; then
    echo "Error: Directory '$track_dir' not found!"
    exit 1
fi
# Start of the XML session
{
    echo '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <Session genome="hg38" locus="ALL" version="8">
        <Resources>'

    # Loop through the bedGraph and bed files in the directory and populate the Resources section
    for file in "$track_dir"/*.bed*; do
        echo "        <Resource path=\"SPLITTED/$(basename $file)\"/>"
    done

    echo '    </Resources>'
    echo '    <Panel height="451" name="DataPanel" width="1783">'

    # Loop through the bedGraph and bed files in the directory and populate the Panel section
    for file in "$track_dir"/*.bed*; do
        base=$(basename "$file" | sed "s/\.bed.*//g")
        # echo "$base"
        if [[ $file == *_1W.bedGraph ]]; then
            echo "        <Track color=\"244,163,97\" strand=\"+1\" id=\"SPLITTED/$(basename $file)\" name=\"$base\" visibility=\"default\" windowFunction=\"mean\">"
            echo "            <DataRange type=\"LINEAR\" absoluteMax=\"100\" absoluteMin=\"0\" relativeMax=\"10\" relativeMin=\"0\"/>"
            echo "        </Track>"
        elif [[ $file == *_2C.bedGraph ]]; then
            echo "        <Track color=\"102,139,138\" strand=\"-1\" id=\"SPLITTED/$(basename $file)\" name=\"$base\" visibility=\"default\" windowFunction=\"mean\">"
            echo "            <DataRange type=\"LINEAR\" absoluteMax=\"100\" absoluteMin=\"0\" relativeMax=\"10\" relativeMin=\"0\"/>"
            echo "        </Track>"
        elif [[ $file == *_3SVstringent.bed ]]; then
            echo "        <Track id=\"SPLITTED/$(basename $file)\" name=\"$base\" displayMode=\"EXPANDED\">"
            echo "            <DataRange maximum=\"100.0\" minimum=\"0.0\" type=\"LINEAR\"/>"
            echo "        </Track>"
        fi
    done

    echo '    </Panel>'
    echo '</Session>'
} >"$output_file"
