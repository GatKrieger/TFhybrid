% plot mutation frequency and changes in binding to all factors

wantedMinScoreTh = -1;
% addition = '_oneOrMoreMuts';
addition = [];
% compareNonFixed = 'sameSite'; 
% compareNonFixed = 'bwSites';

% load(['mat_files\absSumPerPosMutAll', num2str(wantedMinScoreTh), addition, '.mat']);
% load(['mat_files\LRperPosMutAll', num2str(wantedMinScoreTh), addition, '.mat']);
dataDir = [homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\binding_site_mutations\'...
    'mat_files_onlyOneMut_inVitro_std', stdStr, '/'];
load([dataDir, 'absSumPerPosMutAll-1.mat']);
load([dataDir, 'LRperPosMutAll-1.mat']);

flankingBases = (length(fields(LRperPosMutAll.REB1)) - 7)/2;
factors2 = factors;
factors2(ismember(factors2, {'RPB9', 'FHL1', 'STB3', 'PHO2', 'HAP4', 'HMS2', 'FKH2', 'TEC1', 'MIG2'})) = [];
% factors2(ismember(factors2, {'RPB9'})) = [];
NTF = length(factors2);

figure; 
tl = tiledlayout(1, 8);
bgColor = [.9 .9 .9];
cmB = cbrewer2('Blues');
fs=10;
%% PFM
maxLen = max(cellfun(@length, struct2cell(PFMtoUse)));
PFMa = cell(NTF, maxLen);
rectPos4 = nan(NTF, 4);
rectPos4f = nan(NTF, 4);
rectPos2 = nan(NTF, 4);
rectPos3 = nan(NTF, 4);
rectPos1 = nan(NTF, 4);
Ncats = 4; err = 0.5;
for i = 1:NTF
    currMat_allOpts = PFMtoUse.(factors2{i});
    motif_length = length(currMat_allOpts);
    gapSize = maxLen - motif_length;
    halfGapSize = floor(gapSize/2);
    currMat2 = nan(4, maxLen);
    currMat2(:, halfGapSize+1 : halfGapSize+motif_length) = currMat_allOpts;
    PFMa{i} = currMat2;
    rectPos4(i, :) = [halfGapSize + err, err + i*4 - 4, motif_length, 4];
    rectPos4f(i, :) = [halfGapSize + flankingBases + err, err + i*4 - 4, motif_length, 4];
    rectPos2(i, :) = [halfGapSize + flankingBases+ err, err + i*2 - 2, motif_length, 2];
    rectPos3(i, :) = [halfGapSize + flankingBases+ err, err + i*3 - 3, motif_length, 3];
    rectPos1(i, :) = [halfGapSize + flankingBases+ err, err + i - 1, motif_length, 1];
end
PFMa = cell2mat(PFMa);

ax1 = nexttile();
i1 = imagesc(PFMa); hold on;
set(i1, 'alphadata', ~isnan(PFMa)); set(gca, 'color', bgColor);
colormap(ax1, cbrewer2('Greys'));
cb = colorbar('southoutside'); cb.Label.String = 'frequency'; cb.FontSize = fs;
title('Position Frequency Matrix');
for i = 1:NTF; rectangle('Position', rectPos4(i, :)); end
ytickPos = 2.5:Ncats:NTF*Ncats;
caxis([0 1]);
set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
%% consensus sequence
consensusSeq = struct;
Wth = 0.25;
for i = 1:NTF
    currTF = factors2{i};
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
    consNumOrig2 =  [nan(1, flankingBases), consNumOrig, nan(1, flankingBases)];
    consensusSeq.consSeq.(currTF) = origConsSeq2;
    consensusSeq.consNum.(currTF) = consNumOrig2;
end
%% loop
mutFreq_allOpts = cell(NTF, 1);
mutFreq_trans = cell(NTF, 1);
transitionBases = [3,4,1,2];
Nalternatives = cell(NTF, 1);
NalternativesPerNuc = cell(NTF, 1);
deltaBinding_mean = cell(NTF, 1);
deltaBinding_std = cell(NTF, 1);
deltaBinding_trans = cell(NTF, 1);
fixedNonFixed = cell(NTF, 1);
maxLen2 = maxLen + flankingBases*2;
pvals = cell(NTF, 1);
for i = 1:NTF
    currTF = factors2{i};
    currF = LRperPosMutAll.(currTF);
    f = fields(currF);
    currMat_allOpts = nan(4, length(f)); 
    currMat_trans = nan(2, length(f)); 
    currNalternatives = nan(1, length(f));
    currdeltaBinding_mean = nan(4, length(f)); 
    currdeltaBinding_std = nan(4, length(f)); 
    currdeltaBinding_trans = nan(2, length(f));
    currFixedNonFixed = nan(2, length(f));
    currNalternativesPerNuc = nan(4, length(f));
    currPvals = nan(2, length(f));
    currMotif = PFMtoUse.(currTF);
    for j = 1:length(f)
        cPos = f{j};
        currPos = LRperPosMutAll.(currTF).(cPos).bwSp;
        NsitesPerNuc = cellfun(@length, currPos);
        NsitesPerNuc(1) = []; % gap
        currCons = consensusSeq.consNum.(currTF)(j);
        if isnan(currCons)
            [~, currCons] = max(NsitesPerNuc);
        end
        NsitesPerNuc(currCons) = nan;
        NsitesPerNucPrc = NsitesPerNuc ./ nansum(NsitesPerNuc) * 100';
        currMat_allOpts(:, j) = NsitesPerNucPrc;
        currNalternativesPerNuc(:, j) = NsitesPerNuc';
        currTransition = transitionBases(currCons);
        currTransversions = 1:4; currTransversions(currCons) = nan; currTransversions(currTransition) = nan;
        currTransversions(isnan(currTransversions)) = [];
        sumTransversions = nansum(NsitesPerNucPrc(currTransversions));
        currMat_trans(:, j) = [NsitesPerNucPrc(currTransition); sumTransversions];
        currNalternatives(j) = nansum(NsitesPerNuc);
        currPos = LRperPosMutAll.(currTF).(cPos).bwSp(2:5);
        currdeltaBinding_mean(:, j) = cellfun(@mean, currPos);
        currdeltaBinding_std(:, j) = cellfun(@std, currPos);
        currMean = currdeltaBinding_mean(:, j);
        currdeltaBinding_trans(:, j) = [currMean(currTransition); mean(vertcat(currPos{currTransversions}))];
        switch compareNonFixed
            case 'sameSite'
                currFixedNonFixed(:, j) = [absSumPerPosMutAll.(currTF).(cPos).inGenConsLR_allAlt;...
                    LRperPosMutAll.(currTF).(cPos).bwSpConsLR_allAlt];
                currPvals(:, j) = [absSumPerPosMutAll.(currTF).(cPos).inGenConsPval_allAlt;...
                    LRperPosMutAll.(currTF).(cPos).bwSpConsPval_allAlt];
            case 'bwSites'
                currFixedNonFixed(:, j) = [absSumPerPosMutAll.(currTF).(cPos).inGenConsLR_allAlt;...
                    absSumPerPosMutAll.(currTF).(cPos).bwSpConsLR_allAlt];
                currPvals(:, j) = [absSumPerPosMutAll.(currTF).(cPos).inGenConsPval_allAlt;...
                    absSumPerPosMutAll.(currTF).(cPos).bwSpConsPval_allAlt];
        end
    end
    currMat_allOpts2 = nan(4, maxLen2); currMat_trans2 = nan(2, maxLen2); currNalternatives2 = nan(1, maxLen2);
    currNalternativesPerNuc2 = nan(4, maxLen2);
    gapSize = maxLen2 - length(currMat_allOpts);
    halfGapSize = floor(gapSize / 2);
    posToFill = halfGapSize+1 : halfGapSize+length(currMat_allOpts);
    currMat_allOpts2(:, posToFill) = currMat_allOpts;
    currMat_trans2(:, posToFill) = currMat_trans;
    currNalternatives2(:, posToFill) = currNalternatives;
    currNalternativesPerNuc2(:, posToFill) = currNalternativesPerNuc;
    mutFreq_allOpts{i} = currMat_allOpts2; 
    mutFreq_trans{i} = currMat_trans2;
    Nalternatives{i} = currNalternatives2;
    NalternativesPerNuc{i} = currNalternativesPerNuc2;
    currdeltaBinding_mean2 = nan(4, maxLen2); currdeltaBinding_std2 = nan(4, maxLen2); currdeltaBinding_trans2 = nan(2, maxLen2);
    currdeltaBinding_mean2(:, posToFill) = currdeltaBinding_mean;
    deltaBinding_mean{i} = currdeltaBinding_mean2;
    currdeltaBinding_std2(:, posToFill) = currdeltaBinding_std;
    deltaBinding_std{i} = currdeltaBinding_std2;
    currdeltaBinding_trans2(:, posToFill) = currdeltaBinding_trans;
    deltaBinding_trans{i} = currdeltaBinding_trans2;
    currFixedNonFixed2 = nan(2, maxLen2);
    currFixedNonFixed2(:, posToFill) = currFixedNonFixed;
    fixedNonFixed{i} = currFixedNonFixed2;
    currPvals2 = nan(2, maxLen2);
    currPvals2(:, posToFill) = currPvals;
    pvals{i} = currPvals2;
end
%% frequency of cons->alternative : all three options
ax2 = nexttile();
m = cell2mat(mutFreq_allOpts);
i2 = imagesc(m); hold on;
set(i2, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
colormap(ax2, cbrewer2('Blues'));
cb = colorbar('southoutside'); cb.Label.String = 'frequency of alternative'; cb.FontSize = fs;
title('frequency of alternative');
for i = 1:NTF; rectangle('Position', rectPos4f(i, :)); end
ytickPos = 2.5:Ncats:NTF*Ncats; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
%% frequency transition / transversion
ax3 = nexttile();
m = cell2mat(mutFreq_trans);
i3 = imagesc(m); hold on;
set(i3, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
colormap(ax3, cbrewer2('Blues'));
cb = colorbar('southoutside'); cb.Label.String = 'frequency of alternative'; cb.FontSize = fs;
title('transition/transversion');
for i = 1:NTF; rectangle('Position', rectPos2(i, :)); end
ytickPos = 1.5:2:NTF*2; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
%% # site with alternative allele
ax4 = nexttile();
m = cell2mat(Nalternatives);
i4 = imagesc(m); hold on;
set(i4, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
cb = colorbar('southoutside'); cb.Label.String = '# sites'; cb.FontSize = fs;
title('# sites with alternative allele');
for i = 1:NTF; rectangle('Position', rectPos1(i, :)); end
set(gca, 'ytick', 1:NTF, 'yticklabels', factors2);
%% change in binding: cons->alternative, all options
ax5 = nexttile();
m = cell2mat(deltaBinding_mean);
i5 = imagesc(m); hold on;
set(i5, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
cb = colorbar('southoutside'); cb.Label.String = 'log2 ratio alt / cons'; cb.FontSize = fs;
title('\Deltabinding, cons -> alternative');
cm = flipud(cbrewer2('BrBg'));
colormap(ax5, cm);
caxis([-5 5]);
for i = 1:NTF; rectangle('Position', rectPos4f(i, :)); end
ytickPos = 2.5:Ncats:NTF*Ncats; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
%% change in binding: cons->alternative, all options, standard deviation
ax6 = nexttile();
m = cell2mat(deltaBinding_std);
i6 = imagesc(m); hold on;
set(i6, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
cb = colorbar('southoutside'); cb.Label.String = {'log2 ratio alt / cons', 'standard dev'}; cb.FontSize = fs;
title('\Deltabinding, cons -> alternative');
for i = 1:NTF; rectangle('Position', rectPos4f(i, :)); end
ytickPos = 2.5:Ncats:NTF*Ncats; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
% %% change in binding: transition / transversion
% ax7 = nexttile();
% m = cell2mat(deltaBinding_trans);
% i7 = imagesc(m); hold on;
% set(i7, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
% cb = colorbar('southoutside'); cb.Label.String = {'log2 ratio alt / cons'};
% title('\Deltabinding, cons -> alternative');
% colormap(ax7, cm);
% caxis([-5 5]);
% for i = 1:NTF; rectangle('Position', rectPos2(i, :)); end
%% fixed / non-fixed
ax7 = nexttile();
m = cell2mat(fixedNonFixed);
i7 = imagesc(m); hold on;
set(i7, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
cb = colorbar('southoutside'); cb.Label.String = 'log2 ratio alt / cons'; cb.FontSize = fs;
title('Fixed vs. non-fixed mutation');
colormap(ax7, cm);
caxis([-5 5]);
for i = 1:NTF; rectangle('Position', rectPos2(i, :)); end
ytickPos = 1.5:2:NTF*2; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
set(gcf, 'color', 'w');
%% fixed / non-fixed p-values
ax8 = nexttile();
m = cell2mat(pvals);
m1 = m(1:2:NTF*2, :); % fixed mutations
m2 = m(2:2:NTF*2, :); % non-fixed mutations
[h1, crit_p, adj_ci_cvrg, adj_p1] = fdr_bh(m1);
[h2, crit_p, adj_ci_cvrg, adj_p2] = fdr_bh(m2);
H12 = nan(2*NTF, maxLen2);
H12(1:2:NTF*2, :) = h1;
H12(2:2:NTF*2, :) = h2;
i8 = imagesc(H12); hold on;
set(i8, 'alphadata', ~isnan(m)); set(gca, 'color', bgColor);
cb = colorbar('southoutside'); cb.Label.String = 'p-value < 0.05'; cb.FontSize = fs; cb.Ticks = [0, 1];
colormap(ax8, [1 1 1; 0 0 0]);
title('Fixed vs. non-fixed mutation');
for i = 1:NTF; rectangle('Position', rectPos2(i, :)); end
ytickPos = 1.5:2:NTF*2; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
set(gcf, 'color', 'w');