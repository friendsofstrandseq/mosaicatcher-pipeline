BEGIN {OFS="\t"}
$2 != c {prev=0}
NR>1 {
        print $2, prev*s+1,($3+1)*s;
        prev=$3+1;
        c=$2
    }
