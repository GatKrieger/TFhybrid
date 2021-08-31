%% inner loop: calculate substitution matrix and binding consequences for a single position

currPosStr = ['pos', num2str(currPos)];

currPosition = [mats{1}(:, currPos), mats{2}(:, currPos)];
currPosition = currPosition + 1;
allPos = currPosition(:);
currSubMat = nan(Nopts);
currLR = nan(Nopts);
currLRarray = cell(Nopts);
LRperCase = cell(1,Nopts);
currSum = cell(1,Nopts);
absSum = cell(Nopts);
absSum2 = cell(Nopts);
whichSite = cell(Nopts);
currConsSeq = consSeq2{currPos};
if isempty(currConsSeq) | currConsSeq == ' '
    [~, consLoc] = max(histcounts(currPosition, 1:6));
else
    consLoc = find(nucsPlusGap == currConsSeq);
end

for j = 1:5
    for k = 1:5
        op1 = find(currPosition(:, 1) == j & currPosition(:, 2) == k);
        op2 = find(currPosition(:, 1) == k & currPosition(:, 2) == j);
        [uniqueIdx, ii] = unique([op1; op2]);
        currSubMat(j, k) = length(uniqueIdx);
        whichSite{j, k} = uniqueIdx;
        LR = [logRatioCertoPar(op1); logRatioPartoCer(op2)];
        LR = LR(ii);
        currLR(j, k) = nanmean(LR);
        currLRarray{j, k} = LR;
        m = [relevantPeaks(op1, 1); relevantPeaks(op2, 2)];
        mK = [relevantPeaks(op2, 1); relevantPeaks(op1, 2)];
        %             m = m(ii); %
        m(m == 0) = 1; mK(mK == 0) = 1;
        absSum{j, k} = m;
        absSum2{j, k} = mK;
    end
    idx = find(allPos == j);
    currSum{j} = allPeaks(idx);
    
    op1 = find(currPosition(:, 1) == consLoc & currPosition(:, 2) == j);
    op2 = find(currPosition(:, 1) == j & currPosition(:, 2) == consLoc);
    [uniqueIdx, ii] = unique([op1; op2]);
    LR = [logRatioCertoPar(op1); logRatioPartoCer(op2)];
    LR = LR(ii);
    LR = -1.*LR; % second to first
    LR(LR > LRhighTh) = LRhighTh;
    LR(LR < -LRhighTh) = -LRhighTh;
    LRperCase{j} = LR;
end
%%
currSubMat2 = currSubMat';
currSubMat2 = triu(currSubMat2);
subMat{currPos} = currSubMat2';

LRMat(:, currPos) = LRperCase';

occNum = tril(currSubMat);
Nsites = sum(sum(occNum));
freqAltInGen = nan(1, Nopts);
freqAltBwSp = nan(1, Nopts);
LRAltInGen = nan(2, Nopts);
LRAltBwSp = nan(2, Nopts);
absSum_altInGen = nan(2, Nopts);
absSum_altBwSp = nan(2, Nopts);
%%
for j = 1:5
    freqAltInGen(j) = currSubMat(j, j) ./ Nsites * 100;
    freqAltBwSp(j) = currSubMat(j, consLoc) ./ Nsites *100;
    whichSitePerPos.(currPosStr).inGen{j} = whichSite{j,j};
    whichSitePerPos.(currPosStr).bwSp{j} = whichSite{j,consLoc};
    
    LRAltInGen(:, j) = [nanmean(currLRarray{j,j}); nanstd(currLRarray{j,j})];
    LRAltBwSp(:, j) = [nanmean(currLRarray{j,consLoc}); nanstd(currLRarray{j,consLoc})];
    LRperPos.(currPosStr).inGen{j} = currLRarray{j,j};
    LRperPos.(currPosStr).bwSp{j} = currLRarray{j,consLoc};
    
    absSum_altInGen(:, j) = [nanmean(absSum{j,j}); nanstd(absSum{j,j})];
    absSum_altBwSp(:, j) = [nanmean(absSum{j,consLoc}); nanstd(absSum{j,consLoc})];
    absSumPerPos.(currPosStr).inGen{j} = absSum{j,j};
    absSumPerPos.(currPosStr).bwSp{j} = absSum{j,consLoc};
    absSumPerPos.(currPosStr).inGenConsLR{j} = mean(absSum{j,j}) - mean(absSum{consLoc, consLoc});
    [h, p, ci, stats] = ttest2(absSum{j,j}, absSum{consLoc, consLoc});
    absSumPerPos.(currPosStr).inGenConsPval{j} = p;
end
%%
altLoc = 1:5; altLoc(consLoc)= [];
allAlts = [];
allAltBwSp = [];
allLRBwSp = [];
for j = 1:length(altLoc)
    allAlts = [allAlts; absSum{altLoc(j), altLoc(j)}];
    allAltBwSp = [allAltBwSp; absSum{altLoc(j), consLoc}];
    allLRBwSp = [allLRBwSp; currLRarray{altLoc(j), consLoc}];
end
absSumPerPos.(currPosStr).inGenConsLR_allAlt = mean(allAlts) - mean(absSum{consLoc, consLoc});
[h, p, ci, stats] = ttest2(allAlts, absSum{consLoc, consLoc});
absSumPerPos.(currPosStr).inGenConsPval_allAlt = p;

absSumPerPos.(currPosStr).bwSpConsLR_allAlt = mean(allAltBwSp) - mean(absSum{consLoc, consLoc});
[h, p, ci, stats] = ttest2(allAltBwSp, absSum{consLoc, consLoc});
absSumPerPos.(currPosStr).bwSpConsPval_allAlt = p;

allCons_inSNPBwSp = allAltBwSp - allLRBwSp;
LRperPos.(currPosStr).bwSpConsLR_allAlt = mean(allLRBwSp);
[h, p, ci, stats] = ttest2(allCons_inSNPBwSp, allAltBwSp);
LRperPos.(currPosStr).bwSpConsPval_allAlt = p;

abs_altInGen(:, currPos) = absSum_altInGen(1, :)';
abs_altBwSp(:, currPos) = absSum_altBwSp(1, :)';
abs_altInGen_std(:, currPos) = absSum_altInGen(2, :)';
abs_altBwSp_std(:, currPos) = absSum_altBwSp(2, :)';