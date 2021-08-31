% effect of fixed mutations in the motif
indelInFlanking = 1;
mutationRegime = 'onlyOne';
switch mutationRegime
    case 'onlyOne'
        onlyOneMut = 1;
        moreThanOneMut = 0;
    case 'MoreThanOne'
        onlyOneMut = 0;
        moreThanOneMut = 1;
    case 'oneOrMore'
        onlyOneMut = 0;
        moreThanOneMut = 0;
end
        
ifUniqueSite = 1;
% MinPeakHeightByStdfromRandom = 1.5;
flankingBases = 5;
Wth = 0.25;
Nopts = 5; % including indels
peaks_vs_motifScore_params;
origPFM = PFMtoUse.(currTF);
motif_length = size(origPFM, 2);
[maxW, mi] = max(origPFM); 
origConsSeq = nucs(mi);
origConsSeq(maxW <= Wth) = 'N';
consNumOrig = mi;
consNumOrig(maxW <= Wth) = nan;
NwantedPositions = flankingBases*2 + motif_length;
idxOfMotif = flankingBases+1:flankingBases+motif_length;
origConsSeq2 = [repmat('N', 1, flankingBases), origConsSeq, repmat('N', 1, flankingBases)];
MinPeakHeightByStdfromRandomStr = strrep(num2str(MinPeakHeightByStdfromRandom), '.', '_');

% load([homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\peak_tables_std', ...
%     MinPeakHeightByStdfromRandomStr, '\peaks_', currTF, '.mat']);

load([homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\peak_tables_std', ...
    MinPeakHeightByStdfromRandomStr, '_inVitro\peaks_', currTF, '.mat']);

Npeaks = size(peaks, 1);
if ~exist('minScoreToUse'); minScoreToUse = minScore; end
pfmMatchScoreTh = minScoreToUse.(currTF);
ifFigure=0;
absSumPerPosMut = struct;
LRperPosMut = struct;
whichSitePerPosMut = struct;
rawData = struct;
peakTables = struct;
%% loop over all positions
for IDX = 1:NwantedPositions
    cPos = ['pos', num2str(IDX)];
    pfm = origPFM;
    varPos = [];
    if find(ismember(idxOfMotif, IDX))
        % position i is variable, can take any letter
        varPos = IDX - flankingBases;
        pfm(:, varPos) = repmat(.25, 4, 1); 
    end
    bs_mutations3;
    bs_mutations_speciesAgnostic3;
    absSumPerPosMut.(cPos) = absSumPerPos.(cPos);
    LRperPosMut.(cPos) = LRperPos.(cPos);
    whichSitePerPosMut.(cPos) = whichSitePerPos.(cPos);
    rawData.(cPos) = table(mats{1}, mats{2}, relevantPeaks(:, 1), relevantPeaks(:, 2), ...
        topPeaksTable.pfmMatchScorePerSp(:, 1), topPeaksTable.pfmMatchScorePerSp(:, 2), topPeaksTable.gene, ...
        topPeaksTable.correctedLocs, ...
        'variableNames', {'seq_cer', 'seq_par', ...
        'log2_sum_cer', 'log2_sum_par', 'matchScore_cer', 'matchScore_par', 'gene', 'correctedLocs'}); 
    peakTables.(cPos) = topPeaksTable;
end
%%
tVals = nan(Nopts, NwantedPositions);
pVals = nan(Nopts, NwantedPositions);
for i = 1:NwantedPositions
    cPos = ['pos', num2str(i)];
    currBwSp = absSumPerPosMut.(cPos).bwSp;
    currInGen = absSumPerPosMut.(cPos).inGen;
    for j = 1:Nopts
        [h, p, ci, stats] = ttest2(currInGen{j}, currBwSp{j});
        tVals(j,i) = stats.tstat;
        pVals(j,i) = p;
    end
end
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pVals);
pVals = adj_p;
