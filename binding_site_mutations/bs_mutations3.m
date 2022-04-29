% fix the problem of gaps in the ends of the sequence!
%% initiate: 
Npeaks = size(peaks, 1);
motifScoreTh = 0;
% ifUniqueSite = 1;

totalBasesAroundSite_forVis = 301;
basesAroundSite_forVis = (totalBasesAroundSite_forVis-1)/2;
% smallDist = 15; % for analysis
smallDist = motif_length;
Wth = 0.25;
minNsites = 20;

L = size(signalAligned.cer, 1);
motif_length = size(pfm, 2);
[maxW, mi] = max(pfm);
nucs = 'ACGT';
nucsPlusGap = '-ACGT';
consSeq = nucs(mi);
consSeq(maxW <= Wth) = 'N';
consNum = mi;
consNum(maxW <= Wth) = nan;
if exist('origConsSeq')
    consSeq = origConsSeq;
    consNum = consNumOrig;
end

topPeaksTable = peaks;
Npeaks = size(topPeaksTable, 1);
%% align peaks by adjacent maximal motif score, plot aligned peaks 
topPeaksTable = sortrows(topPeaksTable, 'locs');
Ntop = length(topPeaksTable.locs);

if size(signalAligned.cer, 2) ~= length(factors)
    load([homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\data_structs_swalign\signalAligned.mat']);
end
currIdx = find(strcmp(factors, currTF));
currSig = signalAligned.cer(:, currIdx)';

topPeaks = nan(Ntop, totalBasesAroundSite_forVis);
seqCer = cell(Ntop, 1);
seqPar = cell(Ntop, 1);
alignSeq = cell(Ntop, 1);
locationInAlign = nan(Ntop, 1);
[~, maxScoreInSp] = max([topPeaksTable.max_PWMscore_cer, topPeaksTable.max_PWMscore_par], [], 2);
% [~, maxScoreInSp] = max([topPeaksTable.max_score_cer, topPeaksTable.max_score_par], [], 2);

for i = 1:Ntop
    currMaxScoreSp = maxScoreInSp(i);
    currMaxScoreLoc = topPeaksTable.(['max_PWMscore_loc_', sp{currMaxScoreSp}])(i);
    currLoc = topPeaksTable.locs(i) - basesAroundPeak + currMaxScoreLoc - 1;
    smallRegion = currLoc - smallDist : currLoc + smallDist-1;
    smallRegion(smallRegion > L) = [];
    locationInAlign(i) = smallRegion(1);
    seqCer{i} = alignedSeq.cer(smallRegion);
    seqPar{i} = alignedSeq.par(smallRegion);
    alignSeq{i} = wholeAlignment(:, smallRegion);
    
    currRegion = currLoc - basesAroundSite_forVis : currLoc + basesAroundSite_forVis;
    if ~isempty(find(currRegion < 1)) | ~isempty(find(currRegion > L)); continue; end;
    topPeaks(i, :) = currSig(currRegion);
end
sequences = {seqCer, seqPar};

xlabels = -basesAroundSite_forVis:50:basesAroundSite_forVis;
xlabels = cellfun(@num2str, num2cell(xlabels), 'uniformoutput', false);

if ifFigure
    f1 = figure;
    spRange = [10, 2];
    spLocs = {1:2:17, 19, 2:2:18, 20};
    subplot(spRange(1), spRange(2), spLocs{1}); h = imagesc(topPeaks);
    set(h, 'alphadata', ~isnan(topPeaks));
    caxis([0 1000]);
    subplot(spRange(1), spRange(2), spLocs{2});
    plot(nanmean(topPeaks)); axis tight;
    set(gca, 'xtick', 1:50:totalBasesAroundSite_forVis, 'xticklabels', xlabels); set(gcf, 'color', 'w')
end
%% align motifs by PWM

% orderedBS = nan(Ntop, 2);
% alignSeq2 = cell(Ntop, 1);
seqLen = motif_length + smallDist*2;
alignLogical = nan(Ntop, seqLen);
seqCer2 = cell(Ntop, 1);
seqPar2 = cell(Ntop, 1);
betterSeq = cell(Ntop, 1);
correctedLocs = nan(Ntop, 1);
alignSeq2n = cell(Ntop, 1);
s = horzcat(sequences{:});
a = 1:Ntop;
inds = [a' maxScoreInSp];
s2 = sub2ind([Ntop 2], a', maxScoreInSp);
currSequences = s(s2);
pfmMatchScore = nan(Ntop, 1);
Sites1 = locationInAlign;

align_seq_by_PWM;

topPeaksTable.seqDistMotif = motif_length - sum(alignLogical(:, smallDist+1 : smallDist + motif_length), 2);
% flankingCoor = [1:smallDist, smallDist+1+motif_length:size(alignLogical, 2)];
flankingCoor = [smallDist - flankingBases + 1: smallDist, ...
    smallDist + 1 + motif_length: smallDist + motif_length + flankingBases];
topPeaksTable.seqDistFlank = length(flankingCoor) - sum(alignLogical(:, flankingCoor), 2);
seqCer2(find(cellfun(@isempty, seqCer2))) = {''};
seqPar2(find(cellfun(@isempty, seqPar2))) = {''};
topPeaksTable.gcCer = cell2mat(cellfun(@(x) length(regexp(x, 'G|C'))./length(x), seqCer2, 'uniformoutput', false));
topPeaksTable.gcPar = cell2mat(cellfun(@(x) length(regexp(x, 'G|C'))./length(x), seqPar2, 'uniformoutput', false));
topPeaksTable.correctedLocs = correctedLocs;
topPeaksTable.pfmMatchScorePerSp = pfmMatchScorePerSp;
topPeaksTable.strand = Strand;

placeInMat = smallDist+1 : smallDist+motif_length;
placeInMat2 = placeInMat(1)- flankingBases :placeInMat(end)+ flankingBases;
seqCer3 = cellfun(@(x) x(placeInMat2), seqCer2, 'uniformoutput', false);
seqPar3 = cellfun(@(x) x(placeInMat2), seqPar2, 'uniformoutput', false);
topPeaksTable.seq_cer = seqCer3;
topPeaksTable.seq_par = seqPar3;
topPeaksTable.alignLogical = alignLogical(:, placeInMat2);
topPeaksTable_orig = topPeaksTable;
%% mutation from consensus: for best motif score per species
% mut_from_cons_best_motifScore_per_species;
%% unique sites
if ifUniqueSite
    [~, uniqueSites] = unique(correctedLocs, 'stable');
    uniqueSites = intersect(uniqueSites, goodSites);
    alignSeq2n = alignSeq2n(uniqueSites);
    topPeaksTable.row_num = [1:Npeaks]';
    topPeaksTable = topPeaksTable(uniqueSites, :);
    alignLogical = alignLogical(uniqueSites, :);
    pfmMatchScore = pfmMatchScore(uniqueSites);
    betterSeq = betterSeq(uniqueSites);
    disp(['# unique sites = ', num2str(length(uniqueSites))]);
    if length(uniqueSites) < minNsites
        error('no valid sequences');
    else
        disp(['# sites = ', num2str(Npeaks)]);
    end
end

%% plot sequence alignment per site (sweater)
if ifFigure
    currAlign = cell2mat(alignSeq2n);
    M = nuc2num2(currAlign);
    colors = [rgb('White'); rgb('Green'); rgb('Cyan'); rgb('Purple'); rgb('Red')];
    s3 = subplot(spRange(1), spRange(2), spLocs{3});
    imagesc(M);
    colormap(s3, colors);
%     set(gca, 'ytick', 1:3:Ntop*3, 'yticklabels', topPeaksTable.gene)
    subplot(spRange(1), spRange(2), spLocs{4});
    plot(nanmean(alignLogical), 'o-'); axis tight;
    hold on; plot([smallDist+1 smallDist+motif_length], [.8 .8], 'k');
%     ylim([0.65 1]);
    set(gcf, 'color', 'w', 'position', [520    42   691   774]);
end
%% check
% check_betterSeq;
%% plot effect on binding

% bs_mutations_plotEffectOnBinding

%% sustitution matrix

currSeqLen = length(placeInMat2);
consSeq2 = [repmat({''},1, flankingBases), num2cell(consSeq), repmat({''},1, flankingBases)];
topPeaksTable.NmutCore = length(placeInMat) - sum(alignLogical(:, placeInMat), 2);
flnking = alignLogical(:, [placeInMat(1)- flankingBases : placeInMat(1)-1 , ...
    placeInMat(end)+1: placeInMat(end)+ flankingBases]);
topPeaksTable.NmutFlanking = flankingBases*2 - sum(flnking, 2);
