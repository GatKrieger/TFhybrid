# TFhybrid
Code related to the Article: Evolution of transcription factor binding through sequence variations and turnover of binding sites

This is a matlab code that was written for the analysis of TF binding varaition between orthologous genomes in a yeast hybrid.

Experiment: ChEC-seq for 27 TFs in a yeast hybrid (*S. cerevisiae x S. paradoxus*)

Input: bedgraphs from GEO (accession number GSE196451)

The script 'read_data_to_dataStructs.m' generates:
1. Useful data structers
2. Motif enrichment analysis
3. Alignment of orthologous promoters
4. Maps of signal at motifs (poteintial binding sites)
5. Find peaks and associate peaks to motifs
6. Analyis of binding cost due to mutations in the motif
7. Promoter classification
