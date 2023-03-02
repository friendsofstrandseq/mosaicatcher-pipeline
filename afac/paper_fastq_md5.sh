for file in LCL RPE1-WT RPE-BM510 C7; do;
    cd "$file"/fastq
    for f in *.gz; do;
    if [ ! -f "$f".md5 ]
    then
        # echo "File not found"
        echo "$f"
        md5sum "$f" > "$f".md5
    fi
    done;
    cd ../..
done;