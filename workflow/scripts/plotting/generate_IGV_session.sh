#!/usr/bin/env bash

track_dir="$1"
output_file="$2"
reference="$3"

# Debug print statements
printf "Debug: Track directory set to '%s'\n" "$track_dir"
printf "Debug: Output file set to '%s'\n" "$output_file"
printf "Debug: Reference genome set to '%s'\n" "$reference"

# Check if the directory exists
if [[ ! -d "$track_dir" ]]; then
    printf "Error: Directory '%s' not found!\n" "$track_dir" >&2
    exit 1
fi

# Start of the XML session
{
    echo '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <Session genome="'$reference'" locus="ALL" version="8">
        <Resources>'

    # Loop through the bedGraph and bed files in the directory and populate the Resources section
    for file in "$track_dir"/*.bed*; do
        printf "Debug: Processing file '%s'\n" "$file"
        echo "        <Resource path=\"SPLITTED/$(basename $file)\"/>"
    done

    echo '    </Resources>'
    echo '    <Panel height="451" name="DataPanel" width="1783">'

    # Loop through the bedGraph and bed files in the directory and populate the Panel section
    for file in "$track_dir"/*.bed*; do
        base=$(basename "$file" | sed "s/\.bed.*//g")
        printf "Debug: Base name for file '%s' is '%s'\n" "$file" "$base"

        if [[ $file == *_1W.bedGraph ]]; then
            printf "Debug: Matching pattern '*_1W.bedGraph' for file '%s'\n" "$file"
            echo "        <Track color=\"244,163,97\" strand=\"+1\" id=\"SPLITTED/$(basename $file)\" name=\"$base\" visibility=\"default\" windowFunction=\"mean\">"
            echo "            <DataRange type=\"LINEAR\" absoluteMax=\"100\" absoluteMin=\"0\" relativeMax=\"10\" relativeMin=\"0\"/>"
            echo "        </Track>"
        elif [[ $file == *_2C.bedGraph ]]; then
            printf "Debug: Matching pattern '*_2C.bedGraph' for file '%s'\n" "$file"
            echo "        <Track color=\"102,139,138\" strand=\"-1\" id=\"SPLITTED/$(basename $file)\" name=\"$base\" visibility=\"default\" windowFunction=\"mean\">"
            echo "            <DataRange type=\"LINEAR\" absoluteMax=\"100\" absoluteMin=\"0\" relativeMax=\"10\" relativeMin=\"0\"/>"
            echo "        </Track>"
        elif [[ $file == *_3SVstringent.bed ]]; then
            printf "Debug: Matching pattern '*_3SVstringent.bed' for file '%s'\n" "$file"
            echo "        <Track id=\"SPLITTED/$(basename $file)\" name=\"$base\" displayMode=\"EXPANDED\">"
            echo "            <DataRange maximum=\"100.0\" minimum=\"0.0\" type=\"LINEAR\"/>"
            echo "        </Track>"
        else
            printf "Debug: No matching pattern for file '%s'\n" "$file"
        fi
    done

    echo '    </Panel>'
    echo '</Session>'
} >"$output_file"

# Final debug print
printf "Debug: XML session written to '%s'\n" "$output_file"
