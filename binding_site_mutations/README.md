## Mutation analysis scripts:
### Script name: effect_of_fixated_mut_ver3
•	Input: TF name\
•	Code:\
o	Loop over all motif positions, including flanking sequence:\
o	Mutate the specific position so the probabilities in this position are all equal 0.25\
o	Run: bs_mutations3\
o	Run: bs_mutations_speciesAgnostic3\
o	Save outputs in structs

### Script name: bs_mutations3
•	Input: \
o	peaks: list of peaks of a specific TF\
o	pfm: PWM of the specific TF\
•	Code:\
o	Align sequences by the position of the best PWM match\
o	Align sequences to the plus strand of the PWM\
o	Plot sequence alignment per site (optional)\
•	Output:\
o	TopPeaksTable: The original peak table + motif and flanking sequence annotation

### Script name: bs_mutations_speciesAgnostic3
•	Input:\
o	TopPeaksTable\
•	Code:\
o	Take sites where a variant (between orthologues) appears only in the variable position\
o	Report on the number of such sites:\
"# sites for substitution analysis: 770"\
o	Loop over all positions and run: bsmut_innerloop_ver3\
•	Output:\
o	absSumPerPos\
o	LRperPos

### Script name: bsmut_innerloop_ver3
•	Input: \
o	Seq_mats: cell of two matrices, one for each orthologue, with sequences represented by numbers\
o	relevantPeaks: matrix of two vectors, log2 sum of signal at each peak for cer and for par\
o	Generate a substitution matrix per position: the cost in binding of each mutation

