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

peak_prcTh = 95;
peak_table_dir = [homeDir, 'checSeq_project\analyze\new_version\peak_tables_fimo_th\', ...
    'peak_tables_prc', num2str(peak_prcTh), '_fimo/'];
load([peak_table_dir, 'peaks_', currTF, '.mat']);

ifUniqueSite = 1;
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

Npeaks = size(peaks, 1);
if ~exist('minScoreToUse'); minScoreToUse = minScore; end
pfmMatchScoreTh = minScoreToUse.(currTF);
calcMatchScorePerSp=1;
ifFigure=0;
absSumPerPosMut = struct;
LRperPosMut = struct;
whichSitePerPosMut = struct;
rawData = struct;
peakTables = struct;
%% loop over all positions
if all_positions
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
        raw_table = table(seq_mats{1}, seq_mats{2}, relevantPeaks(:, 1), relevantPeaks(:, 2), ...
            topPeaksTable.pfmMatchScorePerSp(:, 1), topPeaksTable.pfmMatchScorePerSp(:, 2), topPeaksTable.gene, ...
            topPeaksTable.correctedLocs, ...
            'variableNames', {'seq_cer', 'seq_par', ...
            'log2_sum_cer', 'log2_sum_par', 'matchScore_cer', 'matchScore_par', 'gene', 'correctedLocs'});
        rawData.(cPos) = raw_table;
        topPeaksTable.seq_cer = seqCer3;
        topPeaksTable.seq_par = seqPar3;
        peakTables.(cPos) = topPeaksTable;
    end
elseif only_pos1
    pfm = origPFM;
    bs_mutations3;
end