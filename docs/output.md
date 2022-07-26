# Outputs

This document describes the final outputs produced by the pipeline. Most of the plots are taken from report generated from the [full-sized test dataset](https://sandbox.zenodo.org/record/1074721) for the pipeline.

The files listed below will be created in the results directory. All paths are relative to the top-level results directory.



## Plots folder

### Mosaic count - reads density across bins

---
File path: `<OUTPUT_FOLDER>/plots/<SAMPLE>/counts/CountComplete.pdf`

Report category: `Mosaic counts`

---



The plots present in this file will allow you to visualize both global statistics at the sample level but also individual Strand-Seq karyotypes at the cell level.

The first page of the pdf file shows the statistics of the analysis result, such
as the distribution of total number of reads per cell, duplication rate, or excluded bins per chromosomes. 


| ![summary](images/plots/stats.png) |
| :--------------------------------: |
| *Global statistics - Sample level* |

Afterwards, every pages show the overview of binning count result of each of the single-cells as presented below. The depth of Crick reads are depicted in the green color in the right side, and the depth of Watson reads are depicted in the orange color in the left side of each chromosome lines. HMM automatically defines the WW/WC/CC status according the reads distribution (yellow background: WC, green background: CC, orange background: WW).


|   ![summary](images/plots/414.png)   |
| :----------------------------------: |
| *Strand-seq karyotype visualisation* |

*Strand-seq karyotype visualisation based on reads counting according defined window (here 100kb). Additionnal statistics are also presented in the upper part of the figure.* 


### SV calls
---
File path: `<OUTPUT_FOLDER>/plots/<SAMPLE>/sv_calls/<CELL>.pdf`

Report category: `SV Calls`


---


### SV consistency
---
File path: `<OUTPUT_FOLDER>/plots/<SAMPLE>/sv_consistency/<CELL>.pdf`

Report category: `SV Consistency`


---


### SV clustering
---
File path: `<OUTPUT_FOLDER>/plots/<SAMPLE>/sv_clustering/<CELL>.pdf`

Report category: `SV Clustering`

---

## Statistics
---
File path: `<OUTPUT_FOLDER>/stats/<SAMPLE>/<CELL>.pdf`

Report category: `Stats`

---

## Raw SV calls (tab-seperated file)
---
File path: `<OUTPUT_FOLDER>/mosaiclassifier/sv_calls/<SAMPLE>/<CELL>.pdf`

---