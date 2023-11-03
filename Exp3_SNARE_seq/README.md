## The SNARE-seq data can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074
| Data type | File name                                                  | Access date  |
| :-------- | :--------------------------------------------------------- | :----------- |
| RNA       | GSE126074_AdBrainCortex_SNAREseq_cDNA.barcodes.tsv.gz      | Nov 3, 2021  |
| RNA       | GSE126074_AdBrainCortex_SNAREseq_cDNA.counts.mtx.gz        | Nov 3, 2021  |
| RNA       | GSE126074_AdBrainCortex_SNAREseq_cDNA.genes.tsv.gz         | Nov 14, 2021 |
| ATAC      | GSE126074_AdBrainCortex_SNAREseq_chromatin.barcodes.tsv.gz | Nov 14, 2021 |
| ATAC      | GSE126074_AdBrainCortex_SNAREseq_chromatin.counts.mtx.gz   | Nov 14, 2021 |
| ATAC      | GSE126074_AdBrainCortex_SNAREseq_chromatin.peaks.tsv.gz    | Nov 14, 2021 |


## We follow the processing step of R package "Signac", who needs download the following externel data： 

| Data name                 | File name                                                                  | Access date  |
| :------------------------ | :------------------------------------------------------------------------- | :----------- |
| fragments.sort.bed.gz     | https://signac-objects.s3.amazonaws.com/snareseq/fragments.sort.bed.gz     | Nov 14, 2021 |
| fragments.sort.bed.gz.tbi | https://signac-objects.s3.amazonaws.com/snareseq/fragments.sort.bed.gz.tbi | May 18, 2018 |
| allen_brain.rds           | https://signac-objects.s3.amazonaws.com/allen_brain.rds                    | Jun 20, 2020 |

And most of the code of S01_Signac_vignette.R was downloaded from https://satijalab.org/signac/articles/snareseq.html.


