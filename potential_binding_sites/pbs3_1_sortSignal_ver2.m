% potential_binding_sites2
% promoter region (400 bp upstream to TSS)

generateNewPeakTables = 0;
fs=10;
peaks_vs_motifScore_params;

PrRegionsCer = find(motifStruct.cer.promotersRegion);
PrRegionsCerAligned = find(ismember(alignedLocV.cer, PrRegionsCer));
PrRegionsPar = find(motifStruct.par.promotersRegion);
PrRegionsParAligned = find(ismember(alignedLocV.par, PrRegionsPar));
PrRegionAligned = {PrRegionsCerAligned, PrRegionsParAligned};
NbestMotifs = 10;
minNsites = 40;
basesAroundPeakForPlot = 40;
basesAroundPeakForNuc = 150;
if ~exist('alMotifScores')
%     score_motif_by_kmers_loop;
    load([homeDir, 'checSeq_project/analyze/signal_on_aligned_seq/data_structs_swalign/alMotifScoresS.mat']);
    load([homeDir, 'checSeq_project/analyze/signal_on_aligned_seq/data_structs_swalign/KmerScoreS.mat']);
end

if size(signalAlignedSm.cer, 2) ~= length(factors)
    load([homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\data_structs_swalign\signalAlignedSm.mat']);
    load([homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\data_structs_swalign\signalAligned.mat']);
end
factors2 = factors;
factors2(ismember(factors2, {'RPB9', 'MIG2', 'TEC1'})) = [];

% reorderedfactors2 = {'ABF1', 'REB1', 'TOD6', 'RAP1', 'TEC1', 'FKH1', 'SUT1', 'DOT6', 'SWI5', 'MIG1', 'MIG2', ...
%     'MSN2', 'SWI4', 'GCR1', 'SKN7', 'SOK2', 'FKH2', 'RGT1', 'ACE2', 'PHO4', 'GCR2', 'STB3', 'PHD1', 'TBF1', ...
%     'HAP4', 'MBP1', 'HMS2', 'FHL1', 'ASH1', 'PHO2'};

NTF = length(factors2);
Sigma = [0, 1, 1.5, 2, 3];
prcVec = 100:-1:1;
prcs = nan(NTF, length(prcVec));

rng(1);
sortedSitesSignal = struct;
sortedSitesIdx = struct;
potentialBSs = struct;
randomSites = struct;
signalRandSites = struct;
sigmasIdx = nan(length(factors2), length(Sigma));
signalAtIdx = nan(length(factors2), length(Sigma));
sigmasPerc = nan(length(factors2), length(Sigma));
sigmasIdx100 = nan(length(factors2), length(Sigma));
sigmasPerc100 = nan(length(factors2), length(Sigma));
for IDX = 1:NTF
    currTF = factors2{IDX};    
    idxTF = find(strcmp(factors, currTF));
    alMotifScores = alMotifScoresS.(currTF);
    KmerScore = KmerScoreS.(currTF);
    currSigSm = [signalAlignedSm.cer(:, idxTF), signalAlignedSm.par(:, idxTF)];
    currSig = [signalAligned.cer(:, idxTF), signalAligned.par(:, idxTF)];
    currCerPBS = [pbsPwm.(currTF).loc_cer_plus_alG; pbsPwm.(currTF).loc_cer_minus_alG]; 
    currParPBS = [pbsPwm.(currTF).loc_par_plus_alG; pbsPwm.(currTF).loc_par_minus_alG]; 
    potentialBS = [currCerPBS; currParPBS];
    whichSp = [ones(length(currCerPBS), 1); 2*ones(length(currParPBS), 1)];
    alignPWM = 0;
    get_signal_at_potential_bs;
    BS = [signalCerPar(whichSp==1, 1); signalCerPar(whichSp==2, 2)];
    potentialBSs.(currTF) = potentialBS;
    % random
    potentialBS = [randsample(PrRegionAligned{1}, length(find(whichSp == 1))), ...
        randsample(PrRegionAligned{2}, length(find(whichSp == 2)))];
    get_signal_at_potential_bs;
    randSites = [signalCerPar(whichSp==1, 1); signalCerPar(whichSp==2, 2)];
    [sortdBS, sortdBSIdx] = sort(BS, 'descend');
    sortedSitesSignal.(currTF) = sortdBS;
    sortedSitesIdx.(currTF) = sortdBSIdx;
    randomSites.(currTF) = potentialBS;
    signalRandSites.(currTF) = randSites;
    sigmaRandSites = std(randSites);
    meanRandSites = mean(randSites);
    
    Points1 = [1, meanRandSites + Sigma(1)*sigmaRandSites, ...
        meanRandSites + Sigma(2)*sigmaRandSites, ...
        meanRandSites + Sigma(3)*sigmaRandSites, ...
        meanRandSites + Sigma(4)*sigmaRandSites];
    Points2 = [Points1(2:end), meanRandSites + Sigma(5)*sigmaRandSites];
    signalAtIdx(IDX, :) = Points2;
    currPrc = prctile(sortdBS, prcVec);
    prcs(IDX, :) = currPrc;
    for sigmas = 1:length(Points2)
        idx = find(sortdBS > Points1(sigmas) & sortdBS <= Points2(sigmas));
        idx = min(idx);
        if isempty(idx); continue; end
        sigmasIdx(IDX, sigmas) = idx;
        sigmasPerc(IDX, sigmas) = idx / (NBS*2) *100;
        
        idx2 = find(currPrc > Points1(sigmas) & currPrc <= Points2(sigmas));
        idx2 = max([0 min(idx2)]);
        sigmasIdx100(IDX, sigmas) = idx2;
        sigmasPerc100(IDX, sigmas) = idx2 / (NBS*2) *100;
    end
end
%% plot sorted signal, + 0-3 standard deviations from random sites
BW = cbrewer2('Greys', 20); BW = BW(end:-1:1, :);
Nsites = structfun(@length, sortedSitesSignal);
[~, si1_5] = sort(sigmasIdx100(:, 3), 'descend');
reorderedfactors2 = factors2(si1_5);

sigmasIdx_2 = [sigmasIdx(:, [3 1]), Nsites];
sigmasIdx_1 = [ones(NTF, 1), sigmasIdx_2(:, 1:end-1)+1];

lw = [1 1 3 1 1];
figure;
tl = tiledlayout(1, 12);

ax1 = nexttile([1 2]);
SigmaToUse = [1 3];
imagesc(prcs(si1_5, :)); hold on;
colormap(ax1, BW); 
% caxis([0 16]);
colors = lines(length(SigmaToUse));
sigmasIdx100_toUse = sigmasIdx100(:, SigmaToUse);
% colors = cbrewer2('OrRd', length(Sigma));
for i = 1:length(factors2)
    Points = sigmasIdx100_toUse(si1_5(i), :);
    yRange = i-.5:i+.5;
    for j = 1:length(SigmaToUse)
        plot([Points(j) Points(j)], yRange, 'color', colors(j, :), 'linewidth', 3);
    end
%     text(100, i, num2str(Nsites(si(i)), '%10.0e'));
end
plot(xlim, [1.5:NTF+.5; 1.5:NTF+.5]', 'w');
set(gca, 'ytick', 1:length(factors2), 'yticklabels', factors2(si1_5));
xlabel('% potential binding sites');
cb = colorbar('southoutside'); cb.Label.String = 'log2 max signal at site'; cb.FontSize = fs;
%% # sites
ax2 = nexttile();
Nsites2 = log10(Nsites(si1_5));
imagesc(Nsites2);
cb = colorbar('southoutside'); cb.Label.String = '# sites'; cb.FontSize = fs;
cb.Ticks = [3 4]; cb.TickLabels = {'10^3', '10^4'};
cm2 = cbrewer2('Blues'); cm2 = cm2(end:-1:1, :);
colormap(ax2, cm2);
set(gca, 'xtick', [], 'ytick', 1:NTF, 'yticklabels', reorderedfactors2);