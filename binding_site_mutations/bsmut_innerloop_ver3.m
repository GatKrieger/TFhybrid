%% inner loop: calculate substitution matrix and binding consequences for a single position

currPosStr = ['pos', num2str(currPos)];

currPosition = [seq_mats{1}(:, currPos), seq_mats{2}(:, currPos)];
currPosition = currPosition + 1;
allPos = currPosition(:);
currSubMat = crosstab(currPosition(:, 1), currPosition(:, 2));
currLR = nan(Nopts);
currLRarray = cell(Nopts);
LRperCase = cell(1,Nopts);
currSum = cell(1,Nopts);
absSum = cell(Nopts);
absSum2 = cell(Nopts);
abs_per_peak = cell(1, Nopts);
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
        if j == k
            op2 = [];
            whichSite{j,j} = op1;
        else
            whichSite{j, k} = op1;
            whichSite{k, j} = op2;
        end
        uniqueIdx = [op1; op2];
        currSubMat(j, k) = length(uniqueIdx);
        LR = [logRatioCertoPar(op1); logRatioPartoCer(op2)];
        currLR(j, k) = nanmean(LR);
        currLRarray{j, k} = LR;
        m = [relevantPeaks(op1, 1); relevantPeaks(op2, 2)];
        mK = [relevantPeaks(op2, 1); relevantPeaks(op1, 2)];
        %             m = m(ii); %
%         m(m == 0) = 1; mK(mK == 0) = 1;
        absSum{j, k} = m;
        absSum2{j, k} = mK;
    end
    idx = find(allPos == j);
    currSum{j} = allPeaks(idx);
    
    cons_in_cer = find(currPosition(:, 1) == consLoc & currPosition(:, 2) == j);
    cons_in_par = find(currPosition(:, 1) == j & currPosition(:, 2) == consLoc);
    if j == consLoc
        % cons_in_cer equals cons_in_par so we need only one.
        cons_in_par = [];
    end
%     whichSite{j} = [cons_in_cer; cons_in_par];
    % LR = log ratio b/w orthologues
    LR = [logRatioCertoPar(cons_in_cer); logRatioPartoCer(cons_in_par)];
    LR = -1.*LR; % alternative to consensus
    LR(LR > LRhighTh) = LRhighTh;
    LR(LR < -LRhighTh) = -LRhighTh;
    LRperCase{j} = LR;
    
    % relevantPeaks: abs binding per orthologue
    curr_abs_cons = [relevantPeaks(cons_in_cer, 1); relevantPeaks(cons_in_par, 2)];
    curr_abs_alt = [relevantPeaks(cons_in_cer, 2); relevantPeaks(cons_in_par, 1)];
    abs_per_peak{j} = [curr_abs_cons, curr_abs_alt];
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
    whichSitePerPos.(currPosStr).inGen{j} = whichSite{j, j};
    whichSitePerPos.(currPosStr).bwSp_alt_cer{j} = whichSite{j, consLoc};
    whichSitePerPos.(currPosStr).bwSp_alt_par{j} = whichSite{consLoc, j};
    
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
[h, p, ci, stats] = ttest2(allAlts, absSum{consLoc, consLoc}, 'tail', 'left');
absSumPerPos.(currPosStr).inGenConsPval_allAlt = p;

absSumPerPos.(currPosStr).bwSpConsLR_allAlt = mean(allAltBwSp) - mean(absSum{consLoc, consLoc});
[h, p, ci, stats] = ttest2(allAltBwSp, absSum{consLoc, consLoc}, 'tail', 'left');
absSumPerPos.(currPosStr).bwSpConsPval_allAlt = p;

% absolute peak level in sites of unique (non-fixed) mutation
absSumPerPos.(currPosStr).abs_cons = {abs_per_peak{1}(:, 1), ...
    abs_per_peak{2}(:, 1), ...
    abs_per_peak{3}(:, 1), ...
    abs_per_peak{4}(:, 1), ...
    abs_per_peak{5}(:, 1)};
absSumPerPos.(currPosStr).abs_alt = {abs_per_peak{1}(:, 2), ...
    abs_per_peak{2}(:, 2), ...
    abs_per_peak{3}(:, 2), ...
    abs_per_peak{4}(:, 2), ...
    abs_per_peak{5}(:, 2)};
% In cases of a unique mutation, test if the conesnsus allele is bound at a
% lower level than cases of a conserved consensus (for reviewer)
levels = cell2mat(absSumPerPos.(currPosStr).abs_cons');
groups = [repmat({'-'}, size(abs_per_peak{1}, 1), 1); ...
    repmat({'A'}, size(abs_per_peak{2}, 1), 1); ...
    repmat({'C'}, size(abs_per_peak{3}, 1), 1); ...
    repmat({'G'}, size(abs_per_peak{4}, 1), 1); ...
    repmat({'T'}, size(abs_per_peak{5}, 1), 1)];
[p, tbl, stats] = kruskalwallis(levels, groups, 'off');
absSumPerPos.(currPosStr).pval_kw_abs_cons = p;

alt_pos = 1:5;
alt_pos(consLoc) = [];
cons_at_alt = cell2mat(absSumPerPos.(currPosStr).abs_cons(alt_pos)');
cons_at_cons = cell2mat(absSumPerPos.(currPosStr).abs_cons(consLoc)');
absSumPerPos.(currPosStr).bwSpConsLR_allCons = ...
    mean(cons_at_alt) - mean(cons_at_cons);
[h, p, ci, stats] = ttest2(cons_at_alt, cons_at_cons, 'tail', 'left');
absSumPerPos.(currPosStr).bwSpConsLR_allCons_pval = p;
% In cases of a unique mutation, test if the alt allele is bound at a
% lower level than cases of a conserved consensus
levels = cell2mat(absSumPerPos.(currPosStr).abs_alt');
groups = [repmat({'-'}, size(abs_per_peak{1}, 1), 1); ...
    repmat({'A'}, size(abs_per_peak{2}, 1), 1); ...
    repmat({'C'}, size(abs_per_peak{3}, 1), 1); ...
    repmat({'G'}, size(abs_per_peak{4}, 1), 1); ...
    repmat({'T'}, size(abs_per_peak{5}, 1), 1)];
[p, tbl, stats] = kruskalwallis(levels, groups, 'off');
absSumPerPos.(currPosStr).pval_kw_abs_alt = p;
% In cases of a common mutation, test if the alt allele is bound at a
% lower level than cases of a conserved consensus
curr_inGen = absSumPerPos.(currPosStr).inGen;
levels = cell2mat(curr_inGen');
groups = [repmat({'-'}, size(curr_inGen{1}, 1), 1); ...
    repmat({'A'}, size(curr_inGen{2}, 1), 1); ...
    repmat({'C'}, size(curr_inGen{3}, 1), 1); ...
    repmat({'G'}, size(curr_inGen{4}, 1), 1); ...
    repmat({'T'}, size(curr_inGen{5}, 1), 1)];
[p, tbl, stats] = kruskalwallis(levels, groups, 'off');
absSumPerPos.(currPosStr).pval_kw_common_alt = p;


allCons_inSNPBwSp = allAltBwSp - allLRBwSp;
LRperPos.(currPosStr).bwSpConsLR_allAlt = mean(allLRBwSp);
[h, p, ci, stats] = ttest2(allCons_inSNPBwSp, allAltBwSp);
LRperPos.(currPosStr).bwSpConsPval_allAlt = p;

abs_altInGen(:, currPos) = absSum_altInGen(1, :)';
abs_altBwSp(:, currPos) = absSum_altBwSp(1, :)';
abs_altInGen_std(:, currPos) = absSum_altInGen(2, :)';
abs_altBwSp_std(:, currPos) = absSum_altBwSp(2, :)';