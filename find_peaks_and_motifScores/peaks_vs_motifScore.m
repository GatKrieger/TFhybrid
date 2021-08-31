%% define parameters
peaks_vs_motifScore_params;
%% get the smoothed signal
factors = fields(strainList);
load PFMhyb;
load pbsPwmHyb;
alMotifScores = alMotifScoresS.(currTF);
% factors = wantedSamples;
strainList2 = struct2array(strainList);
dataDir = [homeDir, 'checSeq_project/analyze/signal_on_aligned_seq/data_structs_swalign/'];
switch averageRepeats
    case 1
        if size(signalAligned.cer, 2) ~= length(factors)
            load([dataDir, 'signalAligned.mat']);
            load([dataDir, 'signalAlignedSm.mat']);
        end
        currIdx = find(strcmp(factors, currSample));
        currStrainList = strainList.(currSample);
        currIdx2 = find(ismember(samples, currStrainList));
        currTF = currSample;
    case 0
        if size(signalAligned.cer, 2) ~= length(strainList2)
            load([dataDir, 'signalAligned_reps.mat']);
            load([dataDir, 'signalAlignedSm_reps.mat']);
        end
        currIdx = find(strcmp(strainList2, currSample));
        currIdx2 = find(ismember(samples, currSample));
        fsl = fields(strainList);
        for i = 1:length(fsl)
            if find(ismember(strainList.(fsl{i}), currSample))
                currTF = fsl{i};
                disp(currTF);
                disp(currSample);
            end
        end        
end
alSig = struct;
for i = 1:2
    currSig = signalAlignedSm.(sp{i})(:, currIdx)';
%     currSig = signalAligned.(sp{i})(:, currIdx)';
    alSig.(sp{i}) = currSig;
end
clear currSig

stats = struct;

MinPeakHeight = [MinPeakHeightFixed MinPeakHeightFixed];

%% open figure
figure;
spRange = [3 4];
spLocs = {1 2 [5, 6], [9 10], 3 4, 7, 8, 11, 12};
suptitle([currTF, ': ', strrep(currSample, '_', ' ')]);
set(gcf, 'color', 'w', 'position', [143         -83        1625         821]);
%% score motifs
score_motif_by_kmers;

get_PWMscore_per_TF;

get_PWMscoreHyb_per_TF;
%% scale signals
Samples = sp;
if ifScale
    switch scaleBy
        case 'sumOnPromoter'
            alignedPromoters = unique(alignedLocGene)';
            alSig_sop = nan(length(alignedPromoters), length(Samples));
            for i = 1:length(alignedPromoters)
                currGeneIdx = find(alignedLocGene == alignedPromoters(i));
                for j = 1:length(Samples)
                    alSig_sop(i,j) = sum(alSig.(Samples{j})(currGeneIdx));
                end
            end
            x = alSig_sop(:, 1);
            y = alSig_sop(:, 2);
%             x = ds_sum_over_promoter.cer(:, currIdx2);
%             y = ds_sum_over_promoter.par(:, currIdx2);
            x(x==0) = nan; y(y==0) = nan;
            if size(x, 2) > 1
                x = nanmean(x, 2);
                y = nanmean(y, 2);
            end
            [new_x, p, p2] = scaleSignal(x, y);
            new_x2 = alSig.cer .* p(1);
            if abs(p2(1)) < 3
                alSig.cer = new_x2;
            end
        case 'fullData'
            x = alSig.cer;
            y = alSig.par;
            [new_x, p, p2] = scaleSignal(x, y);
            if abs(p2(1)) < 3
                alSig.cer = new_x;
            end
    end
    subplot(spRange(1), spRange(2), spLocs{1});
    scatter(x,y);
    x1 = [0 max(max(x), max(y))];
    hold on; plot(x1, polyval(p, x1));
    xl = xlim; yl = ylim;
    maxLim = max(xl(2), yl(2));
    minLim = min(xl(1), yl(1));
    axis([minLim maxLim minLim maxLim]); axis square;
    text(minLim, maxLim, ['y = ', num2str(round(p(1), 2)), 'x + ', num2str(round(p(2), 2))]);
    xlabel('cer'); ylabel('par');
    title('before scaling');
    
    subplot(spRange(1), spRange(2), spLocs{2});
    scatter(new_x, y);
    hold on; plot(x1, polyval(p2, x1));
    title('after scaling');
    xl = xlim; yl = ylim;
    maxLim = max(xl(2), yl(2));
    minLim = min(xl(1), yl(1));
    axis([minLim maxLim minLim maxLim]); axis square;
    text(minLim, maxLim, ['y = ', num2str(round(p2(1), 2)), 'x + ', num2str(round(p2(2), 2))]);
    xlabel('cer'); ylabel('par');
end
clear x y new_x
%% # reads
currNreads = NreadsTable(strcmp(NreadsTable.Row, currSample), :);
currStr = {['cer: ', num2str(currNreads.cer_not_mito, '%10.1e'), ' reads'], ...
        ['par: ', num2str(currNreads.par_not_mito, '%10.1e'), ' reads']};
    dim = [0.0409 0.8973 0.0927 0.0457];
    t = annotation('textbox', dim , 'String', currStr, 'FitBoxToText', 'on');
%% find peaks, sum score around peak, intersect peaks
peaks_vs_motifScore_general;

%% filter peaks by distance from TSS

peaks_dist_from_TSS;

%% plot motif score vs. signal, each species separately

for i = 1:2
    currSp = sp{i};
    subplot(spRange(1), spRange(2), spLocs{4+i});
    x = peaks.(['max_score_', currSp]);
    y = peaks.(['sum_', currSp]);
    y = log2(y+1);
    idx = find(~isnan(x) & ~isnan(y));
    dscatter(x(idx), y(idx));
    r = corr(x, y, 'rows', 'pairwise');
    xlabel('max motif score'); ylabel('log2 ChEC signal'); 
    title({currSp, ['R = ', num2str(round(r,2))]});
    axis square;
end

diffScore = peaks.sum_score_cer - peaks.sum_score_par;
scaleddiffScore = (peaks.sum_score_cer - peaks.sum_score_par)./max([peaks.sum_score_cer,peaks.sum_score_par], [], 2);
diffSumSignal = peaks.sum_cer - peaks.sum_par;
scaleddiffSumSignal = (peaks.sum_cer - peaks.sum_par)./max([peaks.sum_cer, peaks.sum_par], [], 2);

sumcer1 = peaks.sum_cer; sumcer1(sumcer1==0) = 1;
sumpar1 = peaks.sum_par; sumpar1(sumpar1==0) = 1;
peaks.logRatioSum = log2(sumcer1./ sumpar1);
%% plot signal around motif
basesAroundMotif = 150;
nBest = 50;
xlabels1 = -basesAroundMotif:50:basesAroundMotif;
xlabels2 = 1:50:basesAroundMotif*2+1;
xlabels = cellfun(@num2str, num2cell(xlabels1), 'UniformOutput', false);
for j = 1:2
    currSp = sp{j};
    if strcmp(scoreBy, 'av')
        [bestMotifScores, bestMotifs] = maxk(KmerScore, nBest*2); 
    else
        [bestMotifScores, bestMotifs] = maxk(KmerScore(:, j), nBest*2);
    end
%     idxMotif = find(ismember(alGenomeXmers.(currSp), bestMotifs));
    idxMotif = find(ismember(alMotifScores.(currSp), bestMotifScores));
%     idxMotif = find(ismember(alPWMscore.(currSp), bestMotifScores));
    Nsites = length(idxMotif);
    signalAroundMotif = nan(Nsites, basesAroundMotif*2+1);
    currSig = signalAligned.(currSp)(:, currIdx)';
    for i = 1:Nsites
        currLoc = idxMotif(i);
        currRegion = currLoc - basesAroundMotif : currLoc + basesAroundMotif;
        if ~isempty(find(currRegion < 1)) | ~isempty(find(currRegion > length(currSig))); continue; end;
        currS = currSig(currRegion);
        signalAroundMotif(i, :) = currS;
    end
    subplot(spRange(1), spRange(2), spLocs{6+j});
    plot(nanmean(signalAroundMotif)); axis square; axis tight;
    set(gca, 'xtick', xlabels2, 'xticklabels', xlabels);
    xlabel('bases from motif');
    ylabel('signal');
    title(currSp);
end
%% plot: change in motif score vs. change in binding

fs=10;
maxMotifScore = max([peaks.max_score_cer, peaks.max_score_par], [], 2);
motifScoreTh = 0;

maxScoreLR = log2((peaks.max_score_cer+1) ./ (peaks.max_score_par+1));
peaks.maxScoreLR = maxScoreLR;
peaks.sumScoreLR = log2(peaks.sum_score_cer ./ peaks.sum_score_par);
peaks.maxPWMScoreLR = log2((peaks.max_PWMscore_cer+1) ./ (peaks.max_PWMscore_par+1));
peaks.sumPWMScoreLR = log2((peaks.sum_PWMscore_cer+1) ./ (peaks.sum_PWMscore_par+1));

x = maxScoreLR;
y = peaks.logRatioSum;
zaxis = max([peaks.max_score_cer, peaks.max_score_par], [], 2);
zaxis = log2(zaxis);
subplot(spRange(1), spRange(2), spLocs{9});
scatter(x, y, [], zaxis, 'filled', 'markerEdgeColor', 'k');
r = corr(x, y, 'rows', 'pairwise');
currCorrPeakScore = r;

title(['R = ', num2str(round(r, 2))]);
xlabel('max motif score (log2 ratio)');
ylabel('sum signal (log2 ratio)');
cb = colorbar; cb.Label.String = 'log2 max motif score';
axis square; set(gca, 'fontsize', fs);

s2 = subplot(spRange(1), spRange(2), spLocs{10});
%% same, MA plot
% x = log2(mean([peaks.sum_cer, peaks.sum_par], 2));
x = .5 .* (log2(peaks.sum_cer) + log2(peaks.sum_par));
y = peaks.logRatioSum;
% y = log2(peaks.sum_cer) - log2(peaks.sum_par);
z = maxScoreLR;
scatter(x, y, [], z, 'filled', 'markerEdgeColor', [.7 .7 .7]);
cm = cbrewer2('BrBG'); cm = cm(end:-1:1, :);
colormap(s2, cm);
xlabel('log2 mean sum signal');
ylabel('sum signal (log2 ratio)');
cb = colorbar; cb.Label.String = 'max motif score (log2 ratio)';
axis square; set(gca, 'fontsize', fs);

[glist, z] = zscore_fun2(log2(peaks.sum_cer), log2(peaks.sum_par), 1, 'ifFigure', 0, 'prcTh', 0, 'nbins', 20);
peaks.zscore_sumLR = z;
%% print statistics

fprintf('Out of %.0f peaks ...\n', length(peaks.locs));
stats.Npeaks = length(peaks.locs);
peaks_found_in_both = find(peaks.peak_cer == 1 & peaks.peak_par == 1);
peaks_found_in_cer = find(peaks.peak_cer == 1 & peaks.peak_par == 0);
peaks_found_in_par = find(peaks.peak_cer == 0 & peaks.peak_par == 1);
stats.peaks_found_in_both = length(peaks_found_in_both);
stats.peaks_found_in_cer = length(peaks_found_in_cer);
stats.peaks_found_in_par = length(peaks_found_in_par);

diffBound = intersect(find(abs(peaks.logRatioSum) > LRth), glist);
diffBoundLogical = zeros(Npeaks, 1); diffBoundLogical(diffBound) = 1;
stats.diff_bound_peaks = length(diffBound);
stats.motif_cer_signal_cer = length(find(peaks.logRatioSum > LRth & peaks.maxScoreLR > LRth & diffBoundLogical));
stats.motif_par_signal_par = length(find(peaks.logRatioSum < -LRth & peaks.maxScoreLR < -LRth & diffBoundLogical));
stats.motif_cer_signal_par = length(find(peaks.logRatioSum > LRth & peaks.maxScoreLR < -LRth & diffBoundLogical));
stats.motif_par_signal_cer = length(find(peaks.logRatioSum < -LRth & peaks.maxScoreLR > LRth & diffBoundLogical));

fprintf('# peaks found in both: %.0f\n', stats.peaks_found_in_both);
fprintf('# peaks found only in cer: %.0f\n', stats.peaks_found_in_cer);
fprintf('# peaks found only in par: %.0f\n', stats.peaks_found_in_par);
fprintf('# peaks differentially bound: %.0f\n', stats.diff_bound_peaks);
fprintf('# better motif in cer, stronger binding in cer: %.0f\n', stats.motif_cer_signal_cer);
fprintf('# better motif in par, stronger binding in par: %.0f\n', stats.motif_par_signal_par);
fprintf('# better motif in cer, stronger binding in par: %.0f\n', stats.motif_cer_signal_par);
fprintf('# better motif in par, stronger binding in cer: %.0f\n', stats.motif_par_signal_cer);
%%
clear currSig smSignal
