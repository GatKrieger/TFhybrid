# TFhybrid
Code related to the Article: Evolution of transcription factor binding through sequence variations and turnover of binding sites.\
This is a matlab code for the analysis of TF binding varaition between orthologous genomes in a yeast hybrid.\
Experiment: ChEC-seq for 27 TFs in a yeast hybrid (*S. cerevisiae x S. paradoxus*)\
Input: bedgraphs from GEO (accession number GSE196451)

<img src="https://user-images.githubusercontent.com/60549750/166002522-85e555da-80bb-4bd4-b6f0-cdc0f3152e38.svg" width="200" height="348">

The script 'runall.m' generates:
1. Useful data structers
2. Motif enrichment analysis
3. Alignment of orthologous promoters
4. Maps of signal at motifs (poteintial binding sites)
5. Find peaks and associate peaks to motifs
6. Analyis of binding cost due to mutations in the motif
7. Promoter classification
