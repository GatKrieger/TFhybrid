%% classification params
decide_using_corr = 0;
norm_to_max=1;
norm_to_total=0;
peakStandardToUse = 'background';
whichScore = '7mer';

%% load data 
currDir = [homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\peak_tables_std1_5_inVitro\'];
load([currDir, 'peaks_', currTF, '.mat']);
if ~exist('MinPeakHeightPerTFTable')
    load([homeDir, 'checSeq_project\analyze\potential_binding_sites\MinPeakHeightPerTFTable_inVitro.mat'])
end
Npeaks = length(peaks.gene);
%% filter motif score

Ktop = 5;
maxMScores_PWM = [peaks.max_PWMscore_cer, peaks.max_PWMscore_par];
get_PWMscore_per_TF;
bestMotifs_PWM = maxk(unique(genPWMscore.cer), Ktop);

Ktop = 20;
maxMScores_PWMsHyb = [peaks.max_PWMsHyb_cer, peaks.max_PWMsHyb_par];
get_PWMscoreHyb_per_TF;
bestMotifs_PWMsHyb = maxk(unique(genPWMscore_hyb.cer), Ktop);

Ktop = 20;
maxMScores_7mer = [peaks.max_score_cer, peaks.max_score_par];
bestMotifs_7mer = maxk(unique(alMotifScoresS.(currTF).cer), Ktop);

switch whichScore
    case 'PWM'
        maxMScores = maxMScores_PWM;
        bestMotifs = bestMotifs_PWM;
    case 'PWMsHyb'
        maxMScores = maxMScores_PWMsHyb;
        bestMotifs = bestMotifs_PWMsHyb;
    case '7mer'
        maxMScores = maxMScores_7mer;
        bestMotifs = bestMotifs_7mer;
end
minMotif = min(bestMotifs);
PWMscore_perfect = [peaks.max_PWMscore_cer, peaks.max_PWMscore_par] == 1;
peaks.hasMotif = maxMScores > minMotif | PWMscore_perfect;
peaks.hasMotif_sum = sum(peaks.hasMotif, 2);
hasMotif = find(any(peaks.hasMotif, 2));
minPeakHeight = MinPeakHeightPerTFTable.mean_plus_1_5_std(ismember(MinPeakHeightPerTFTable.factor, currTF));
peakHeight = [peaks.height_cer, peaks.height_par];
peaks.hasPeak = peakHeight > minPeakHeight;
%% find the motif score in the exact location of the maximal score of the other species

% motifSites = nan(Npeaks, 2);
motifScore_bySite_7mer = struct;
motifScore_bySite_7mer.cer = nan(Npeaks, 2);
motifScore_bySite_7mer.par = nan(Npeaks, 2);

motifScore_bySite_PWM = motifScore_bySite_7mer;
motifScore_bySite_PWMsHyb = motifScore_bySite_7mer;

otherSp = {'par', 'cer'};

for i = 1:Npeaks
    for j = 1:2
        currLoc = peaks.locs(i) - basesAroundPeak -1 + peaks.(['max_score_loc_', sp{j}])(i);  
        currScore = alMotifScoresS.(currTF).(sp{j})(currLoc);
        currScore_otherSp = alMotifScoresS.(currTF).(otherSp{j})(currLoc);
        motifScore_bySite_7mer.(sp{j})(i, :) = [currScore currScore_otherSp];
        
        currLoc = peaks.locs(i) - basesAroundPeak -1 + peaks.(['max_PWMscore_loc_', sp{j}])(i);  
        currScore = alPWMscore.(sp{j})(currLoc);
        currScore_otherSp = alPWMscore.(otherSp{j})(currLoc);
        motifScore_bySite_PWM.(sp{j})(i, :) = [currScore currScore_otherSp];
        
        currLoc = peaks.locs(i) - basesAroundPeak -1 + peaks.(['max_PWMsHyb_loc_', sp{j}])(i);  
        currScore = alPWMscore_hyb.(sp{j})(currLoc);
        currScore_otherSp = alPWMscore_hyb.(otherSp{j})(currLoc);
        motifScore_bySite_PWMsHyb.(sp{j})(i, :) = [currScore currScore_otherSp];
    end
end
%% classify peaks by motif score
% ms_bin = motif score binary

ms7mer_bin = [motifScore_bySite_7mer.cer > min(bestMotifs_7mer), ...
    motifScore_bySite_7mer.par > min(bestMotifs_7mer)];
msPWM_bin = [motifScore_bySite_PWM.cer == 1, ...
    motifScore_bySite_PWM.par == 1];

same_PWMloc = peaks.max_PWMscore_loc_cer == peaks.max_PWMscore_loc_par;
same_PWMsHybloc = peaks.max_PWMsHyb_loc_cer == peaks.max_PWMsHyb_loc_par;

peaks.consMotif = (all(ms7mer_bin == 1, 2) | all(msPWM_bin == 1, 2));

peaks.closeTO = ((((ms7mer_bin(:, 1) == 1 & ms7mer_bin(:, 2) == 0 & ...
    ms7mer_bin(:, 3) == 1 & ms7mer_bin(:, 4) == 0) | ...
    (ms7mer_bin(:, 1) == 1 & ms7mer_bin(:, 2) == 1 & ms7mer_bin(:, 3) == 1 & ms7mer_bin(:, 4) == 0) | ...
    (ms7mer_bin(:, 1) == 1 & ms7mer_bin(:, 3) == 1 & ms7mer_bin(:, 4) == 1 & ms7mer_bin(:, 2) == 0)) ...
    & same_PWMsHybloc == 0));

peaks.AddMotif = ((ms7mer_bin(:, 1) == 1 & ms7mer_bin(:, 2) == 0 & ms7mer_bin(:, 3) == 0) | ...
    (ms7mer_bin(:, 3) == 1 & ms7mer_bin(:, 1) == 0 & ms7mer_bin(:, 4) == 0) | ...
    (msPWM_bin(:, 1) == 1 & msPWM_bin(:, 2) == 0 & msPWM_bin(:, 3) == 0) | ...
    (msPWM_bin(:, 3) == 1 & msPWM_bin(:, 1) == 0 & msPWM_bin(:, 4) == 0) );

peaks.closeTO(peaks.consMotif==1) = 0;
peaks.AddMotif(peaks.consMotif==1) = 0;
motifDef = cell(Npeaks, 1);
motifDef(peaks.consMotif==1) = {'consMotif'};
motifDef(peaks.closeTO==1) = {'closeTO'};
motifDef(peaks.AddMotif==1) = {'AddMotif'};
motifDef(peaks.consMotif==0 & peaks.closeTO==0 & peaks.AddMotif==0) = {'noMotif'};
peaks.motifDef = motifDef;
%% classify binding site by peak height & motif score

BSdef = cell(Npeaks, 1);
bothAboveBg = sum(peaks.isAboveBg, 2) == 2;
bothHavePeak = sum(peaks.hasPeak, 2) == 2;
if strcmp(peakStandardToUse, 'background')
    peakStandard = bothAboveBg;
elseif strcmp(peakStandardToUse, 'peak')
    peakStandard = bothHavePeak;
end

hasMotif_cer = peaks.hasMotif(:, 1) == 1;
hasMotif_par = peaks.hasMotif(:, 2) == 1;

consBS = find(peakStandard & peaks.consMotif==1);
TOBS = find(peakStandard & peaks.closeTO==1);
addBS_cer = ((peaks.AddMotif==1 & peaks.isAboveBg(:, 1) == 1 & hasMotif_cer == 1) | ...
    (peaks.closeTO == 1 & peaks.isAboveBg(:, 1) == 1 & hasMotif_cer == 1 & ...
    peaks.isAboveBg(:, 2) == 0));
addBS_par = ((peaks.AddMotif==1 & peaks.isAboveBg(:, 2) == 1 & hasMotif_par == 1) | ...
    (peaks.closeTO == 1 & peaks.isAboveBg(:, 2) == 1 & hasMotif_par == 1 & ...
    peaks.isAboveBg(:, 1) == 0));

addBS = find(addBS_cer | addBS_par);

BSdef(:) = {'noDef'};
BSdef(consBS) = {'cons'};
BSdef(TOBS) = {'TO'};
BSdef(addBS) = {'Addition'};
peaks.BSdef = BSdef;
peaks.addBS_cer = addBS_cer;
peaks.addBS_par = addBS_par;
sums = [peaks.sum_cer, peaks.sum_par];
[max_sum, mi] = max(sums, [], 2);
peaks.max_sum = max_sum;
%% classify promoters

uniqPrs = unique(peaks.gene(hasMotif), 'stable');
prT = cell(length(uniqPrs), 2);
prT(:, 1) = uniqPrs;
for i = 1:length(uniqPrs)
    currPr = uniqPrs{i};
    idxPr = find(ismember(peaks.gene, currPr));
    typeD = peaks.BSdef(idxPr);   
    stype = unique(typeD);
    consBS = 0; TOBS = 0; addBS = 0;
    if find(ismember(stype, 'cons')); consBS = 1; end
    if find(ismember(stype, 'Addition')); addBS = 1; end
    if find(ismember(stype, 'TO')); TOBS = 1; end
    
    hasMotifD = peaks.hasMotif(idxPr, :);
    sHasMotif = sum(hasMotifD, 1);
    
    sumAddBS = sum([peaks.addBS_cer(idxPr), peaks.addBS_par(idxPr)], 1);
    sumAddBS = sumAddBS > 0;
    
    if consBS & ~TOBS & ~addBS
        prT{i, 2} = 'cons';
        
    elseif TOBS | (addBS & sum(sumAddBS) == 2) 
        prT{i, 2} = 'TO';
        
    elseif addBS & consBS 
        prT{i, 2} = 'add';
        
    elseif addBS & ~consBS 
        prT{i, 2} = 'addOnly';
    end
end
prT(cellfun(@isempty, prT(:, 2)), 2) = {'noDef'};
prT = table(prT(:, 1), prT(:, 2), 'variableNames', {'pr', 'type'});
t = tabulate(prT.type)
groups = {'cons', 'TO', 'add', 'addOnly', 'noDef'};

%% correlation of signal and motif score along promoter

peak_turnover;
% decide for turnover using correlation:
% promoters with low signal correlation and low motif correlation are
% turning over.
% promoters with peak addtion only (add only): also filter by low signal correlation.
corrTh = 0.8;
idx = ismember_smart(allGenes, prT.pr);
prT.corr_signal = corrSignalPerPr(idx);
prT.corr_motif = corrMotifPerPr(idx);
diffPattern = prT.corr_signal < corrTh & prT.corr_motif < corrTh;
prT.des = prT.type;
tabulate(prT.des);
%% mark if close or far TO
Npr = length(prT.pr);
motifDefs = cell(Npr, 1);
TOtype = cell(Npr, 1);
BSdef = cell(Npr,1);
Nbs = nan(Npr, 1);
maxPeakSum = cell(Npr, 1);
for i = 1:Npr
    currPr = prT.pr{i};
    idx_peaks = find(strcmp(peaks.gene, currPr));
    
    currBSdef = peaks.BSdef(idx_peaks);
    idxToTake = find(~ismember(currBSdef, {'noDef'}));
    currBSdef = currBSdef(idxToTake);
    BSdef{i} = currBSdef;
    Nbs(i) = length(currBSdef);
    
    currPeaks = peaks(idx_peaks(idxToTake), :);
    maxPeakSum{i} = currPeaks.max_sum;
    
    currmotifDef = currPeaks.motifDef;
    motifDefs{i} = {currmotifDef{:}};
    NcloseTO = find(ismember(currmotifDef, 'closeTO'));
    currAddBS = [peaks.addBS_cer(idx_peaks), peaks.addBS_par(idx_peaks)];
    hasTwoAddedSite = find(sum(currAddBS, 2) == 2);
    if strcmp(prT.des{i}, 'TO')
        if any(NcloseTO)
            TOtype{i} = 'close';
%             if length(NcloseTO) == 2 & range(currPeaks) < totalbasesAroundPeak
%                 BSdef{i} = {'closeTO'};
            s = sum([currPeaks.peak_cer, currPeaks.peak_par]);
            m = max(s);
            Nbs(i) = m;
        elseif any(hasTwoAddedSite > 0)
            TOtype{i} = 'close';
        else
            TOtype{i} = 'far';
            s = sum([currPeaks.peak_cer, currPeaks.peak_par]);
            m = max(s);
            Nbs(i) = m;
        end
    end
end
prT.motifDef = motifDefs;
prT.BSdef = BSdef;
prT.NBS = Nbs;
prT.maxPeakSum = maxPeakSum;
prT.TOtype = TOtype;
%%  analyze, plot
%toPlot = 'corr_signal';
%dToPlot = cell(1, length(groups));
figure; 
spRange = [1 5];
Ng = length(groups);
subplot(spRange(1), spRange(2), 1); hold on;
for i = 1:Ng
    currGidx = find(ismember(prT.des, groups{i}));
    scatter(prT.corr_signal(currGidx), prT.corr_motif(currGidx));
end
% legend(groups, 'location', 'southoutside', 'orientation', 'horizontal');
axis square;
xlabel('corr signal'); ylabel('corr motif');
l = min([xlim, ylim]);
axis([l 1 l 1]);

% sum over promoter:
idxTF = find(ismember(factors, currTF));
sop = [ds_av.cer.sum_over_promoter(:, idxTF), ds_av.par.sum_over_promoter(:, idxTF)];
sop_max = max(sop, [], 2);
% sop_max = log2(sop_max+1);
gidx = ismember_smart(allGenes, prT.pr);
prT.sop_max = sop_max(gidx);
pr_groups = cell(1, length(groups));
for i = 1:Ng
    currGidx = find(ismember(prT.des, groups{i}));
    pr_groups{i} = prT.sop_max(currGidx);
end
subplot(spRange(1), spRange(2), 2);
plotSpread(pr_groups);
set(gca, 'xtick', 1:length(groups), 'xticklabels', groups);
set(gca, 'yscale', 'log');
ylabel('sum signal over promoter (max allele)');
axis square;

subplot(spRange(1), spRange(2), 3); hold on;
peakTypes = {'consMotif', 'closeTO', 'AddMotif'};
peakTypeN = sum([peaks.consMotif, peaks.closeTO, peaks.AddMotif])';
bar(peakTypeN); 
set(gca, 'xtick', 1:length(peakTypes), 'xticklabels', peakTypes);
axis square;
ylabel('# peaks');

subplot(spRange(1), spRange(2), 4); hold on;
prTypeN = nan(1, length(groups));
for i = 1:Ng
    prTypeN(i) = length(find(ismember(prT.des, groups{i})));
end
bar(prTypeN);
set(gca, 'xtick', 1:length(groups), 'xticklabels', groups);
axis square;
ylabel('# promoters');

signalBothAlleles = sum(sop, 2);
% signalBothAlleles = log2(signalBothAlleles+1);
prT.sop_both = signalBothAlleles(gidx);
sum_group = nan(1, length(groups));
for i = 1:Ng
    currGidx = find(ismember(prT.des, groups{i}));
    if norm_to_total
        sum_group(i) = sum(prT.sop_both(currGidx));
    end
    if norm_to_max
        sum_group(i) = sum(prT.sop_max(currGidx));
    end
end
subplot(spRange(1), spRange(2), 5);
if norm_to_total
    prcOfTotalSignal = sum_group ./ sum(prT.sop_both)*100;
end
if norm_to_max
    prcOfTotalSignal = sum_group ./ sum(prT.sop_max)*100;
end
bar(prcOfTotalSignal);
set(gca, 'xtick', 1:length(groups), 'xticklabels', groups);
axis square;
ylabel('% of signal (linear)');

suptitle(currTF);
set(gcf, 'position', [1          41        1600         783], 'color', 'w');

tabulate(prT.TOtype(strcmp(prT.des, 'TO')))
%% peak level per group
maxPeakSum = cell(1, Ng);
groupID = {};
for i = 1:Ng
    currGroup = groups{i};
    idx = find(strcmp(prT.type, currGroup));
    d = vertcat(prT.maxPeakSum{idx});
%     d = log2(d+1);
    maxPeakSum{i} = d;
    groupID = vertcat(groupID, repmat({currGroup}, length(d), 1));
end
medianPeakLevel = cellfun(@median, maxPeakSum);
y = vertcat(maxPeakSum{:});
peakLevelTable = table(y, groupID, 'variableNames', {'maxSum', 'BSdef'});
% figure; plotSpread(maxPeakSum);
% 
% [~, ~, statsAnova] = anova1(y, groupID);
% [c,~,~,gnames] = multcompare(statsAnova);
%%
savedir = 'classified_peaks/';
save([savedir, 'peaks_', currTF, '.mat'], 'peaks');
save([savedir, 'prT_', currTF, '.mat'], 'prT');