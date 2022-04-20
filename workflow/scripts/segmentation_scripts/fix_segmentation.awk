# FIX SEGMENTATION
BEGIN {OFS="\t"} 
    {
        if ($1 == name) { 
            $12 = int(($14-1)/int(window)); 
            }
        print
    }
