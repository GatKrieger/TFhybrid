# TFhybrid
Code related to the Article: Evolution of transcription factor binding through sequence variations and turnover of binding sites.\
This is a matlab code for the analysis of TF binding varaition between orthologous genomes in a yeast hybrid.\
Experiment: ChEC-seq for 27 TFs in a yeast hybrid (*S. cerevisiae x S. paradoxus*)\
Input: bedgraphs from GEO (accession number GSE196451)

![TFhybrid_scheme](https://user-images.githubusercontent.com/60549750/165992154-a56484d6-dce1-49df-9260-1c21566aa34f.jpg)

The script 'runall.m' generates:
1. Useful data structers
2. Motif enrichment analysis
3. Alignment of orthologous promoters
4. Maps of signal at motifs (poteintial binding sites)
5. Find peaks and associate peaks to motifs
6. Analyis of binding cost due to mutations in the motif
7. Promoter classification
