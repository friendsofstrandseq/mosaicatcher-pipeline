#!/usr/bin/awk -f

BEGIN {OFS="\t"}

NR==FNR {
    norm[$1,$2,$3]=$5;
    next
}

FNR==1 {
    if (NF < 10) {
        print $0, "bin_id";
        add_bin_id = 1;
    } else {
        print $0;
    }
}

FNR>1 {
    $8=(($1,$2,$3) in norm && norm[$1,$2,$3]!="good") ? "None" : $8;
    if (add_bin_id) {
        print $0, $1"_"$2"_"$3;
    } else {
        print $0;
    }
}
