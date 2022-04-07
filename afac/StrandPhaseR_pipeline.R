#!/usr/bin/Rscript

args=commandArgs(TRUE)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[6]) )

library(StrandPhaseR)

strandPhaseR(
    inputfolder="/g/korbel2/weber/MosaiCatcher_files/bam_HJ/TALL03-DEA5/selected", 
    outputfolder="/g/korbel2/weber/MosaiCatcher_output_sample_HJ/strand_states/TALL03-DEA5/100000.selected_j0.1_s0.5_scedist20/StrandPhaseR_analysis.chr11", 
    configfile = "/g/korbel2/weber/MosaiCatcher_output_sample_HJ/strand_states/TALL03-DEA5/100000.selected_j0.1_s0.5_scedist20/StrandPhaseR.chr11.config",  
    WCregions = "/g/korbel2/weber/MosaiCatcher_output_sample_HJ/strand_states/TALL03-DEA5/100000.selected_j0.1_s0.5_scedist20/strandphaser_input.txt", 
    positions="/g/korbel2/weber/MosaiCatcher_files/snv_sites_to_genotype/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz", 
    )
