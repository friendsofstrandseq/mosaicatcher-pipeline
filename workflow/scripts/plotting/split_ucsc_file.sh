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
/^track/ {
    if(filename != "") close(filename);

    # Extract track name and modify it for filename
    if($0 ~ /name="([^"]+)"/) {
        filename = gensub(/.*name="([^"]+)".*/, "\\1", 1);
    } else if ($0 ~ /name=([^ ]+)/) {
        filename = gensub(/.*name=([^ ]+).*/, "\\1", 1);
    }

    filename = (filename ~ /_W$/ ? gensub(/_W$/, "_1W.bedGraph", 1, filename) :
                (filename ~ /_C$/ ? gensub(/_C$/, "_2C.bedGraph", 1, filename) :
                (filename ~ /_SV_stringent$/ ? gensub(/_SV_stringent$/, "_3SVstringent.bed", 1, filename) : filename)));

    # Redirect to the specified output directory
    filename = outdir "/" filename;
}
{
    print $0 > filename;
}'

echo "Tracks have been saved to $output_dir/"
