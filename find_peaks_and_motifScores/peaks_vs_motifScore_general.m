%% peaks vs motif score: general (compare species and compare replcates)


%% find peaks: general
Samples = fields(alSig);
peakInfo = struct;
disp(['Min peak height: ', num2str(round(MinPeakHeight, 2))]);
proms = unique(alignedLocGene, 'stable');
ax = cell(2,1);
for i = 1:length(Samples)
    currSample = Samples{i};
    currSig = alSig.(currSample);
    subplot(spRange(1), spRange(2), spLocs{2+i});
    currMinPeakHight = MinPeakHeight(i);
    currMinPeakProminence = currMinPeakHight;
    findpeaks(currSig, 'MinPeakHeight', currMinPeakHight, ...
        'MinPeakDistance', MinPeakDistance,'MinPeakProminence', currMinPeakProminence, ...
        'MinPeakWidth', MinPeakWidth, 'Annotate', 'extents');
    legend off;
    yl = ylim; ylim([0 yl(2)]);
    [pks,locs,w,p] = findpeaks(currSig, 'MinPeakHeight', currMinPeakHight, ....
        'MinPeakDistance', MinPeakDistance, 'MinPeakProminence', currMinPeakProminence, ...
        'MinPeakWidth', MinPeakWidth, 'Annotate', 'extents');
    peakInfo.(currSample) = table(locs', w', round(p'), pks', ...
        'variableNames', {'locs', 'width', 'prominence', 'height'});
    if restrictKpeaks
        currKpeaks = min(Kpeaks, length(p));
        [~, idx] = maxk(p', currKpeaks);
        peakInfo.(currSample) = peakInfo.(currSample)(idx, :);
    end
    lowPeaks = zeros(length(pks), 1);
    for ii = 1:length(proms)
        currRegion = find(alignedLocGene == proms(ii));
        locsInRegion = find(ismember(locs, currRegion));
        if ~isempty(locsInRegion)
            sortdLocs = sort(locs(locsInRegion));
            currRegion = sortdLocs(1) - confRange : sortdLocs(end) + confRange;
            currRegion(currRegion < 1) = []; currRegion(currRegion > L) = [];
            currSigRegion = currSig(currRegion);
            currSigRegion(currSigRegion == 0) = nan;
            goodSigHere = prctile(currSigRegion, prctileToFlt);
            peakDiff = pks(locsInRegion) - goodSigHere;
            idxLow = find(pks(locsInRegion) < goodSigHere);
            lowPeaks(locsInRegion(idxLow)) = 1;
        end
    end
    peakInfo.(currSample)(find(lowPeaks), :) = [];
    fprintf('found %.0f peaks\n', length(peakInfo.(currSample).locs));
    title([currSample, sprintf(': found %.0f peaks\n', length(peakInfo.(currSample).locs))]);
    xlabel('position on aligned promoters');
    ylabel('smoothed chec-seq signal');
end

stats.(['peaks_', Samples{1}]) = length(peakInfo.(Samples{1}).locs);
stats.(['peaks_', Samples{2}]) = length(peakInfo.(Samples{2}).locs);
peakInfo.(Samples{1}).(['peak_', Samples{1}]) = ones(length(peakInfo.(Samples{1}).locs), 1);
peakInfo.(Samples{1}).(['peak_', Samples{2}])= zeros(length(peakInfo.(Samples{1}).locs), 1);
peakInfo.(Samples{2}).(['peak_', Samples{1}]) = zeros(length(peakInfo.(Samples{2}).locs), 1);
peakInfo.(Samples{2}).(['peak_', Samples{2}]) = ones(length(peakInfo.(Samples{2}).locs), 1);
%% signal and motif score around peak

peaks = [peakInfo.(Samples{1}); peakInfo.(Samples{2})];
Npeaks = size(peaks, 1);
signalAroundPeak = nan(Npeaks,2);
sumMScoreAroundPeak = nan(Npeaks,2);
maxMScoreAroundPeak = nan(Npeaks, 2);
MScoreAroundPeakAtLocOfOther = nan(Npeaks, 2);
maxMScoreLocationAroundPeak = nan(Npeaks, 2);
sumPWMScoreAroundPeak = nan(Npeaks,2);
maxPWMScoreAroundPeak = nan(Npeaks, 2);
PWMScoreAroundPeakAtLocOfOther = nan(Npeaks, 2);
maxPWMScoreLocationAroundPeak = nan(Npeaks, 2);
sumPWMScoreHybAroundPeak = nan(Npeaks,2);
maxPWMScoreHybAroundPeak = nan(Npeaks, 2);
PWMScoreHybAroundPeakAtLocOfOther = nan(Npeaks, 2);
maxPWMScoreHybLocationAroundPeak = nan(Npeaks, 2);
height = nan(Npeaks, 2);
otherSample = fliplr(Samples');
for i = 1:Npeaks
    currLoc = peaks.locs(i);
    currRegionLims = [currLoc - basesAroundPeak, currLoc + basesAroundPeak];
    currRegionLims(1) = max([currRegionLims(1), 1]);
    currRegionLims(2) = min([currRegionLims(2), L]);
    currRegionForMotif = currRegionLims(1):currRegionLims(2);
    
    currRegionLimsForSumChecSignal = [currLoc - basesAroundPeakForSum, currLoc + basesAroundPeakForSum];
    currRegionLimsForSumChecSignal(1) = max([currRegionLimsForSumChecSignal(1), 1]);
    currRegionLimsForSumChecSignal(2) = min([currRegionLimsForSumChecSignal(2), L]);
    currRegionForSum = currRegionLimsForSumChecSignal(1):currRegionLimsForSumChecSignal(2);
    
    for j = 1:length(Samples)
        currSample = Samples{j};
        currSignal = alSig.(currSample)(currRegionForSum);
        signalAroundPeak(i, j) = sum(currSignal);
        height(i, j) = alSig.(currSample)(currLoc);
        
        currMSScore = alMotifScores.(currSample)(currRegionForMotif);
        sumMScoreAroundPeak(i, j) = nansum(currMSScore);
        [m, mi] = max(currMSScore);
        maxMScoreAroundPeak(i, j) = m;
        maxMScoreLocationAroundPeak(i, j) = mi;
        MScoreAroundPeakAtLocOfOther(i, j) = alMotifScores.(otherSample{j})(currRegionForMotif(mi));
        
        currPWMscore = alPWMscore.(currSample)(currRegionForMotif);
        sumPWMScoreAroundPeak(i, j) = nansum(currPWMscore);
        [m, mi] = max(currPWMscore);
        maxPWMScoreAroundPeak(i, j) = m;
        maxPWMScoreLocationAroundPeak(i, j) = mi;
        
        currPWMscore = alPWMscore_hyb.(currSample)(currRegionForMotif);
        sumPWMScoreHybAroundPeak(i, j) = nansum(currPWMscore);
        [m, mi] = max(currPWMscore);
        maxPWMScoreHybAroundPeak(i, j) = m;
        maxPWMScoreHybLocationAroundPeak(i, j) = mi;
    end
end
peaks.gene = allGenes(alignedLocGene(peaks.locs));
peaks.geneid = alignedLocGene(peaks.locs)';

peaks.(['height_', Samples{1}]) = height(:, 1);
peaks.(['height_', Samples{2}]) = height(:, 2);
peaks.(['sum_', Samples{1}]) = signalAroundPeak(:, 1);
peaks.(['sum_', Samples{2}]) = signalAroundPeak(:, 2);
peaks.sum_mean = nanmean(signalAroundPeak, 2);
signalAroundPeak(signalAroundPeak == 0) = 1;
peaks.(['log2_sum_', Samples{1}]) = log2(signalAroundPeak(:, 1));
peaks.(['log2_sum_', Samples{2}]) = log2(signalAroundPeak(:, 2));
peaks.log2_sum_mean = nanmean(log2(signalAroundPeak), 2);
peaks.sum_mean = nanmean(signalAroundPeak, 2);

peaks.(['sum_score_', Samples{1}]) = sumMScoreAroundPeak(:, 1);
peaks.(['sum_score_', Samples{2}]) = sumMScoreAroundPeak(:, 2);
peaks.sum_score_mean = nanmean(sumMScoreAroundPeak, 2);
peaks.(['max_score_', Samples{1}]) = maxMScoreAroundPeak(:, 1);
peaks.(['max_score_', Samples{2}]) = maxMScoreAroundPeak(:, 2);
peaks.(['max_score_loc_', Samples{1}]) = maxMScoreLocationAroundPeak(:, 1);
peaks.(['max_score_loc_', Samples{2}]) = maxMScoreLocationAroundPeak(:, 2);

peaks.(['sum_PWMscore_', Samples{1}]) = sumPWMScoreAroundPeak(:, 1);
peaks.(['sum_PWMscore_', Samples{2}]) = sumPWMScoreAroundPeak(:, 2);
peaks.sum_PWMscore_mean = nanmean(sumPWMScoreAroundPeak, 2);
peaks.(['max_PWMscore_', Samples{1}]) = maxPWMScoreAroundPeak(:, 1);
peaks.(['max_PWMscore_', Samples{2}]) = maxPWMScoreAroundPeak(:, 2);
peaks.(['max_PWMscore_loc_', Samples{1}]) = maxPWMScoreLocationAroundPeak(:, 1);
peaks.(['max_PWMscore_loc_', Samples{2}]) = maxPWMScoreLocationAroundPeak(:, 2);

peaks.(['sum_PWMsHyb_', Samples{1}]) = sumPWMScoreHybAroundPeak(:, 1);
peaks.(['sum_PWMsHyb_', Samples{2}]) = sumPWMScoreHybAroundPeak(:, 2);
peaks.sum_PWMsHyb_mean = nanmean(sumPWMScoreHybAroundPeak, 2);
peaks.(['max_PWMsHyb_', Samples{1}]) = maxPWMScoreHybAroundPeak(:, 1);
peaks.(['max_PWMsHyb_', Samples{2}]) = maxPWMScoreHybAroundPeak(:, 2);
peaks.(['max_PWMsHyb_loc_', Samples{1}]) = maxPWMScoreHybLocationAroundPeak(:, 1);
peaks.(['max_PWMsHyb_loc_', Samples{2}]) = maxPWMScoreHybLocationAroundPeak(:, 2);

wholeAlSeq = [alPromSeq.Alignment];
[peaks, siLocs] = sortrows(peaks, 'locs');
%% intersect peaks
% [C, ia, ic] = unique(peaks.locs);
% peaks = peaks(ia, :);
peaks.width = round(peaks.width);
peaks.height = round(peaks.height);
Npeaks = length(peaks.locs);

halfWidth = 10;
locs1 = peaks.locs - halfWidth;
locs2 = peaks.locs + halfWidth;
ifIntersect = zeros(Npeaks, 1);
idxIntersect = zeros(Npeaks, 2);
peakLogical = [peaks.(['peak_', Samples{1}]), peaks.(['peak_', Samples{2}])];
for i = 2:Npeaks
    if locs1(i) <= locs2(i-1)
        idxIntersect(i, :) = [i, i-1];
        peakmax = nan(1, 2);
        peakmax(1) = max(peaks.(['height_', Samples{1}])(i-1), peaks.(['height_', Samples{2}])(i-1));
        peakmax(2) = max(peaks.(['height_', Samples{1}])(i), peaks.(['height_', Samples{2}])(i));
        peakLogical(i, :) = sum([peakLogical(i-1, :); peakLogical(i, :)]);
        if peakmax(1) > peakmax(2)
            ifIntersect(i) = 1;
        else
            ifIntersect(i-1) = 1;
        end
    end
end
peakLogical(peakLogical > 1) = 1;
peaks.(['peak_', Samples{1}]) = peakLogical(:, 1);
peaks.(['peak_', Samples{2}]) = peakLogical(:, 2);
stats.N_intersect_peaks = length(find(ifIntersect));
stats.N_unique_peaks = length(peaks.locs) - stats.N_intersect_peaks;
fprintf('intersecting peaks: %.0f\n', stats.N_intersect_peaks);
peaks(find(ifIntersect), :) = [];

Npeaks = length(peaks.locs);
fprintf('# peaks: %.0f\n', Npeaks);
