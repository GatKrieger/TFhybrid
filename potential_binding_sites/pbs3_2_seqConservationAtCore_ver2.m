%% seq. conservation around core motif
calcMatchScorePerSp = 0;
smallDist = 10; %
nucs = 'ACGT';
BSlocCorrected = struct;
err=0.5;
Nsigmas = size(sigmasIdx_2, 2);
seqConsAtBS = struct;
% minScore2 = minScore; minScore2.TEC1 = -2;
for IDX = 1:NTF
    currTF = reorderedfactors2{IDX};
    currCerPBS = [pbsPwm.(currTF).loc_cer_plus_alG; pbsPwm.(currTF).loc_cer_minus_alG]; 
    currParPBS = [pbsPwm.(currTF).loc_par_plus_alG; pbsPwm.(currTF).loc_par_minus_alG]; 
    potentialBS = [currCerPBS; currParPBS];
    pbsStrand = [ones(length(pbsPwm.(currTF).loc_cer_plus_alG), 1); ...
        zeros(length(pbsPwm.(currTF).loc_cer_minus_alG), 1); ...
        ones(length(pbsPwm.(currTF).loc_par_plus_alG), 1); ...
        zeros(length(pbsPwm.(currTF).loc_par_minus_alG), 1)];
    NBS = length(potentialBS);
    get_PWMscore_per_TF;
    halfMotifLen = floor(motif_length/2);
    seqLen = halfMotifLen*2 + smallDist*2 + 1;
    alignLogical = zeros(NBS, seqLen);
    for i = 1:NBS
        currLoc = potentialBS(i);
        currSite = currLoc - halfMotifLen - smallDist : currLoc + halfMotifLen + smallDist;
        currSite(currSite > L) = []; currSite(currSite<1) = 1;
        currSeq = wholeAlignment(:, currSite);
        currSimilarity = currSeq(2, :);
        currSeq(2, :) = [];
        if pbsStrand(i) == 0
            currSeq(1, :) = seqrcomplement(currSeq(1, :));
            currSeq(2, :) = seqrcomplement(currSeq(2, :));
            currSimilarity = fliplr(currSimilarity);
        end
        tmp = zeros(1, seqLen);
        tmp(currSimilarity == '|') = 1;
        alignLogical(i, :) = tmp;
    end

    al = alignLogical(sortedSitesIdx.(currTF), :);
    m = nan(Nsigmas, size(al, 2));
    currsigmasIdx_1 = sigmasIdx_1(si1_5(IDX), :);
    currsigmasIdx_2 = sigmasIdx_2(si1_5(IDX), :);
    currsigmasIdx_1(isnan(currsigmasIdx_1)) = [];
    currsigmasIdx_2(isnan(currsigmasIdx_2)) = [];
    for J = 1:length(currsigmasIdx_1)
        currIdx = currsigmasIdx_1(J) : currsigmasIdx_2(J);
        if length(currIdx) < 5; continue; end
        currAL = al(currIdx, :);
        m(J, :) = mean(currAL, 1);
    end
    seqConsAtBS.(currTF) = m;
end

%% 
slen = structfun(@length, seqConsAtBS);
maxlen = max(slen);
minlen = min(slen);
Ncats = 3;
seqConsMat = cell(NTF, 1);
rectPos = nan(NTF, 4);

for i = 1:NTF
    currTF = reorderedfactors2{i};
    if ~isfield(seqConsAtBS, currTF) 
        seqConsMat{i} = nan(Ncats, maxlen); 
        continue; 
    end
    currMat = seqConsAtBS.(currTF);
    currMat2 = nan(Ncats, maxlen);
    currMatLen = size(currMat, 2);
    sizeDiff = maxlen - currMatLen;
    gapSize = round(sizeDiff/2);
    currMat2(:, gapSize+1 : gapSize + currMatLen) = currMat;
    seqConsMat{i} = currMat2;
    currMotifLength = length(PFMtoUse.(currTF));
    sizeDiff2 = maxlen - currMotifLength;
    gapSize2 = round(sizeDiff2/2);
    rectPos(i, :) = [gapSize2+.5, .5+i*3-3, currMotifLength, 3];
end
seqConsMat2 = cell2mat(seqConsMat);
% seqConsMat2 = seqConsMat2(:, 6:39);
% rectPos(:, 1) = rectPos(:, 1) - 5;
%% plot coservation motif score
ax3 = nexttile([1 2]);
i3 = imagesc(seqConsMat2); hold on;
caxis([0.6 1]);
set(i3, 'alphadata', ~isnan(seqConsMat2));
hold on;
for i = 1:NTF
    currTF = reorderedfactors2{i};
    if ~isfield(seqConsAtBS, currTF); continue; end;
    rectangle('Position', rectPos(i, :));
end
cb = colorbar('southoutside'); cb.Label.String = 'sequence conservation'; cb.FontSize = fs;
set(gca, 'ytick', Ncats-1:Ncats:NTF*Ncats, 'yticklabels', reorderedfactors2);
lineGaps = Ncats+.5: Ncats : NTF*Ncats+.5;
plot(xlim, [lineGaps; lineGaps]', 'k');
currSz = size(seqConsMat2, 2);
currMid = currSz/2;
currDistances = -15:5:15;
currDistancesStr = cellfun(@num2str, num2cell(currDistances), 'uniformoutput', false);
xVec = currMid + currDistances;
set(gca, 'xtick', xVec, 'xticklabels', currDistancesStr);
xlabel('Distance from motif center (bp)');