if ~exist('tss_wholeGenIdx'); load('tss_wholeGenIdx.mat'); end;
dist_from_TSS = nan(Npeaks, 1);
closestGeneID = nan(Npeaks, 1);
genesWTss = find(~isnan(tss_struct.cer(:, 1)) & ~isnan(tss_struct.par(:, 1)));
tss_cer_wholeGenIdx1 = tss_wholeGenIdx.cer(genesWTss);
allGenes1 = allGenes(genesWTss);

for i = 1:Npeaks
    currPeakLoc_aligned = peaks.locs(i);
    currRealLoc = alignedLocV.cer(currPeakLoc_aligned);
    n=1;
    while isnan(currRealLoc)
        currRealLoc = alignedLocV.cer(currPeakLoc_aligned + n); n = n+1;
    end
    currDiffVec = abs(currRealLoc - n - tss_cer_wholeGenIdx1);
    [m, mi] = min(currDiffVec);
    dist_from_TSS(i) = m;
    closestGeneID(i) = genesWTss(mi);
end

a = horzcat(peaks.geneid, closestGeneID);
[~, idxDiff] = setdiff(peaks.geneid, closestGeneID);

peaks.distFromClosestTSS = dist_from_TSS;
peaks.closestGene = closestGeneID;
uniqueGenes = unique(peaks.gene, 'stable');
for i = 1:length(uniqueGenes)
    idxGene = find(strcmp(peaks.gene, uniqueGenes{i}));
    dist_TSS_currPr = peaks.distFromClosestTSS(idxGene);
    [~, mi] = min(dist_TSS_currPr);
    currClosestGene = closestGeneID(idxGene(mi));
    peaks.closestGene(idxGene) = (currClosestGene);
end

%% check

% gIds = nan(Npeaks, 1);
% for i = 1:Npeaks; gIds(i) = find(ismember(allGenes, peaks.closestGene{i})); end
locsReal = nan(Npeaks, 2);
TSSReal = nan(Npeaks, 2);
for i = 1:Npeaks
    currPeakLoc_aligned = peaks.locs(i);
    for j = 1:2
        currSp = sp{j};
        currRealLoc = alignedLocV.(currSp)(currPeakLoc_aligned);
        n=1;
        while isnan(currRealLoc)
            currRealLoc = alignedLocV.(currSp)(currPeakLoc_aligned + n); n = n+1;
        end
        currRealTSS = tss_wholeGenIdx.(currSp)(peaks.geneid(i));
        locsReal(i, j) = currRealLoc;
        TSSReal(i, j) = currRealTSS;
    end
end

distFromTSS = locsReal - TSSReal -1;
distFromTSS = abs(distFromTSS);
peaks.distFromTSS_cer = distFromTSS(:, 1);
peaks.distFromTSS_par = distFromTSS(:, 2);
peaks.locsReal_cer = locsReal(:, 1);
peaks.locsReal_par = locsReal(:, 2);
peaks.TSSReal_cer = TSSReal(:, 1);
peaks.TSSReal_par = TSSReal(:, 2);
%%
peaks(peaks.distFromClosestTSS > distFromTSSTh, :) = [];
Npeaks = length(peaks.gene);
disp(['# peaks = ', num2str(Npeaks)]);
