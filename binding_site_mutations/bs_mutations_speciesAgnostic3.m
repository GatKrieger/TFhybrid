% binding site mutations: mutation per site, regardless if it came from cer
% or from par

consSeq1_2 = consSeq;
consSeq1_2(strfind(consSeq, 'N')) = ' ';
consSeq2 = [repmat({''},1, flankingBases), num2cell(consSeq1_2), repmat({''},1, flankingBases)];

seqCer3 = cellfun(@(x) x(placeInMat2), seqCer2(uniqueSites), 'uniformoutput', false);
seqPar3 = cellfun(@(x) x(placeInMat2), seqPar2(uniqueSites), 'uniformoutput', false);
seqCer_onlyCore = cellfun(@(x) x(placeInMat), seqCer2(uniqueSites), 'uniformoutput', false);
seqPar_onlyCore = cellfun(@(x) x(placeInMat), seqPar2(uniqueSites), 'uniformoutput', false);
if indelInFlanking == 0
    % do not allow indels at flanking region
    indelCases = find(contains(seqCer3, '-') | contains(seqPar3, '-'));
else
    % Allow indels at flanking region, not in core motif
    indelCases = find(contains(seqCer_onlyCore, '-') | contains(seqPar_onlyCore, '-'));
end
%%
currAlignLogical = topPeaksTable.alignLogical(:, idxOfMotif);
% No variation other than in the variable position
posStrict = 1:length(placeInMat);
posStrict(varPos) = [];
% Allow variability in variable places in the motif (like in abf1) addition 16/01/2022
posStrict(isnan(consNumOrig)) = [];
sameNucInStrictPos = sum(currAlignLogical(:, posStrict), 2);
Nmuts = length(posStrict) - sameNucInStrictPos;
mutAtStrictPos = find(Nmuts > 0);
NoMutAtStrictPos = find(Nmuts == 0);

forSub = 1:length(uniqueSites);

if onlyOneMut & ~moreThanOneMut
    forSub([indelCases; mutAtStrictPos]) = [];
elseif ~onlyOneMut & moreThanOneMut
    forSub([indelCases; NoMutAtStrictPos]) = [];
elseif ~onlyOneMut & ~moreThanOneMut
    forSub(indelCases) = [];
end
disp(['# sites for substitution analysis: ', num2str(length(forSub))]);
topPeaksTable = topPeaksTable(forSub, :);
seqCer3 = seqCer3(forSub);
seqPar3 = seqPar3(forSub);

seqA = {seqCer3, seqPar3};
seq_mats = cell(1, 2);
for i = 1:length(seqA)
    currSeq = cell2mat(seqA{i});
    seq_mats{i} = nuc2num2(currSeq);
end
s = horzcat(seqCer3, seqPar3);
alignSeq2n = alignSeq2n(forSub);

allSubstitutions = cell(4,4);
for i = 1:4
    for j = 1:4
        allSubstitutions{i,j} = [nucs(j), nucs(i)];
    end
end
A = ones(4);
idx = find(triu(A, 1));
allSubstitutions(idx) = {''};
%% mutation effect: both parts

relevantPeaks = [topPeaksTable.sum_cer, topPeaksTable.sum_par];
relevantPeaks(relevantPeaks == 0) = 1; relevantPeaks = log2(relevantPeaks);
allPeaks = relevantPeaks(:);
logRatioCertoPar = topPeaksTable.logRatioSum;
logRatioPartoCer = -1 .* topPeaksTable.logRatioSum;
LRhighTh = 5;

subMat = cell(1, length(placeInMat2));
signalMat = cell(Nopts, length(placeInMat2));
LRMat = cell(Nopts, length(placeInMat2));

LRperPos = struct;
absSumPerPos = struct;
whichSitePerPos = struct;

abs_altInGen = nan(Nopts, length(placeInMat2));
abs_altBwSp = nan(Nopts, length(placeInMat2));
abs_altInGen_std = nan(Nopts, length(placeInMat2));
abs_altBwSp_std = nan(Nopts, length(placeInMat2));
%% loop over all positions
% 
% for i = 1:length(placeInMat2)
%     currPos = i;
%     bs_mutations_speciesAgnostic3_innerLoop;
% end
%% Take only the wanted position

currPos = IDX;
% bs_mutations_speciesAgnostic3_innerLoop;
bsmut_innerloop_ver3;