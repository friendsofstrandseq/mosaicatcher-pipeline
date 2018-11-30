BEGIN { srand(seed); 
        OFS="\t"
    } 
    {
        vaf=min_vaf+rand()*(max_vaf-min_vaf);
        print $0, vaf
    }
