#!/usr/bin/env bash

input_gz_file="$1"

# Check if the file exists
if [ ! -f "$input_gz_file" ]; then
    echo "Error: File '$input_gz_file' not found!"
    exit 1
fi

output_dir="$2"
mkdir -p "$output_dir"

# Process the gzipped file and split it into separate tracks
zcat "$input_gz_file" | awk -v outdir="$output_dir" '
BEGIN { filename = "" }
/^track/ {
    if (filename != "") close(filename);

    # Extract track name and modify it for filename
    if ($0 ~ /name="([^"]+)"/) {
        filename = gensub(/.*name="([^"]+)".*/, "\\1", "g", $0);
    } else if ($0 ~ /name=([^ ]+)/) {
        filename = gensub(/.*name=([^ ]+).*/, "\\1", "g", $0);
    }

    # Apply substitutions based on the suffix
    if (filename ~ /_W$/) {
        filename = gensub(/_W$/, "_1W.bedGraph", "g", filename);
    } else if (filename ~ /_C$/) {
        filename = gensub(/_C$/, "_2C.bedGraph", "g", filename);
    } else if (filename ~ /_SV_stringent$/) {
        filename = gensub(/_SV_stringent$/, "_3SVstringent.bed", "g", filename);
    }

    # Redirect to the specified output directory
    filename = outdir "/" filename;
}
{
    print > filename;
}'

echo "Tracks have been saved to $output_dir/"
