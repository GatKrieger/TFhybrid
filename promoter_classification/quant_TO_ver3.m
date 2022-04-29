%% classification params
decide_using_corr = 0;
norm_to_max = 1;
norm_to_total = 0;
peakStandardToUse = 'background';

%% load data 

peak_tables_file = [homeDir, 'checSeq_project/analyze/new_version/peak_tables_ed.xlsx'];
peaks_mat_dir = [homeDir, 'checSeq_project/analyze/new_version/peak_tables_fimo_th/peak_tables_prc95_fimo/'];
peaks_ed = readtable(peak_tables_file, 'Sheet', currTF, 'PreserveVariableNames', true);
load([peaks_mat_dir, 'peaks_', currTF, '.mat']);

cols_needed = {'peak_cer', 'peak_par', 'motif_cer', 'motif_par', ...
    'motif_com_alt', 'motif_uni_alt', 'motif_com_cons', ...
    'conservedMotifs', 'nonConservedMotifs', 'hasMotif', 'noMotifInBoth'};
peaks(:, cols_needed) = peaks_ed(:, cols_needed);

if flt_peaks
    %     prc_top_peaks = 75;
    peak_height = [peaks.height_cer, peaks.height_par];
    peak_th = prctile(peak_height(:), prc_top_peaks);
    peaks_bin = peak_height >= peak_th;
    below_th = find(all(peaks_bin == 0, 2));
    peaks(below_th, :) = [];
end
Npeaks = height(peaks);

%% genome-wide PWM score

get_PWMscore_per_TF;
get_PWMscoreHyb_per_TF;

%% find the motif score in the exact location of the maximal score of the other species

motifScore_bySite_7mer = struct;
motifScore_bySite_7mer.cer = nan(Npeaks, 2);
motifScore_bySite_7mer.par = nan(Npeaks, 2);

motifScore_bySite_PWM = motifScore_bySite_7mer;
motifScore_bySite_PWMsHyb = motifScore_bySite_7mer;

otherSp = {'par', 'cer'};
basesAroundPeak = 30;

loc_correct = -1;
if strcmp(currTF, 'ABF1')
    loc_correct = 0;
end

for i = 1:Npeaks
    for j = 1:2
        currLoc = peaks.locs(i) - basesAroundPeak + loc_correct + peaks.(['max_score_loc_', sp{j}])(i);  
        currScore = alMotifScoresS.(currTF).(sp{j})(currLoc);
        currScore_otherSp = alMotifScoresS.(currTF).(otherSp{j})(currLoc);
        motifScore_bySite_7mer.(sp{j})(i, :) = [currScore currScore_otherSp];
        
        currLoc = peaks.locs(i) - basesAroundPeak + loc_correct + peaks.(['max_PWMscore_loc_', sp{j}])(i);  
        currScore = alPWMscore.(sp{j})(currLoc);
        currScore_otherSp = alPWMscore.(otherSp{j})(currLoc);
        motifScore_bySite_PWM.(sp{j})(i, :) = [currScore currScore_otherSp];
        
        currLoc = peaks.locs(i) - basesAroundPeak + loc_correct + peaks.(['max_PWMsHyb_loc_', sp{j}])(i);  
        currScore = alPWMscore_hyb.(sp{j})(currLoc);
        currScore_otherSp = alPWMscore_hyb.(otherSp{j})(currLoc);
        motifScore_bySite_PWMsHyb.(sp{j})(i, :) = [currScore currScore_otherSp];
    end
end
%% classify peaks by motif score

% use FIMO definitions
if strcmp(currTF, 'ABF1')
    mscores = [peaks.max_PWMscore_cer, peaks.max_PWMscore_par];
    mscores = unique(mscores);
    mscores(mscores == 0) = [];
    min_score = min(mscores);
else
    curr_ftable = ftable(ismember(ftable.motif_id, currTF), :);
    motif_seq = unique(curr_ftable.matched_sequence);
    pfm = PFMtoUse.(currTF);
    motif_length = size(pfm, 2);
    max_score = prod(max(pfm));
    scores = nan(length(motif_seq), 1);
    for i = 1:length(motif_seq)
        curr_seq = motif_seq{i};
        curr_motif_num = nuc2num2(curr_seq);
        scores(i) = prod(pfm(sub2ind(size(pfm), curr_motif_num, 1:motif_length)));
    end
    scores(scores == 0) = [];
    scores = scores ./ max_score;
    min_score = min(scores);
end
pwm_bin = [motifScore_bySite_PWM.cer >= min_score, ...
    motifScore_bySite_PWM.par >= min_score];
% Use the 7-mer score as well
Ktop = 20;
m7_scores = [peaks.max_score_cer, peaks.max_score_par];
m7_min_score = min(maxk(unique(m7_scores), Ktop));
m7_bin = [motifScore_bySite_7mer.cer >= m7_min_score, ...
    motifScore_bySite_7mer.par >= m7_min_score];
if ismember(currTF, 'REB1')
    m_bin = pwm_bin & m7_bin;
else
    m_bin = pwm_bin;
end

same_PWMloc = peaks.max_PWMscore_loc_cer == peaks.max_PWMscore_loc_par;
same_PWMsHybloc = peaks.max_PWMsHyb_loc_cer == peaks.max_PWMsHyb_loc_par;

peaks.consMotif = (all(m7_bin == 1, 2) | all(pwm_bin == 1, 2));

closeTO_op1 = (...
    m_bin(:, 1) == 1 & m_bin(:, 2) == 0 & ...
    m_bin(:, 3) == 1 & m_bin(:, 4) == 0);
closeTO_op2 = (...
    m_bin(:, 1) == 1 & m_bin(:, 2) == 1 & ...
    m_bin(:, 3) == 1 & m_bin(:, 4) == 0);
closeTO_op3 = (...
    m_bin(:, 1) == 1 & m_bin(:, 2) == 0 & ...
    m_bin(:, 3) == 1 & m_bin(:, 4) == 1);
closeTO_must = same_PWMsHybloc == 0;
closeTO_opts = [closeTO_op1, closeTO_op2, closeTO_op3, closeTO_must];
peaks.closeTO = (closeTO_op1 | closeTO_op2 | closeTO_op3) & closeTO_must;

addMotif_op1 = (...
    m_bin(:, 1) == 1 & m_bin(:, 2) == 0 & ...
    m_bin(:, 3) == 0);
addMotif_op2 = (...
    m_bin(:, 1) == 0 & ...
    m_bin(:, 3) == 1 & m_bin(:, 4) == 0);
addMotif_op3 = (...
    m_bin(:, 1) == 1 & m_bin(:, 2) == 0 & ...
    m_bin(:, 3) == 0);
addMotif_op4 = (...
    m_bin(:, 1) == 0 & ...
    m_bin(:, 3) == 1 & m_bin(:, 4) == 0);
peaks.AddMotif = addMotif_op1 | addMotif_op2 | addMotif_op3 | addMotif_op4;

peaks.closeTO(peaks.consMotif==1) = 0;
peaks.AddMotif(peaks.consMotif==1) = 0;
motifDef = cell(Npeaks, 1);
motifDef(peaks.consMotif==1) = {'consMotif'};
motifDef(peaks.closeTO==1) = {'closeTO'};
motifDef(peaks.AddMotif==1) = {'AddMotif'};
motifDef(peaks.consMotif==0 & peaks.closeTO==0 & peaks.AddMotif==0) = {'noMotif'};
peaks.motifDef = motifDef;
fprintf('Conserved motif: %.0f\nclose TO motifs: %.0f\nAdd motifs: %.0f\nNo motif: %.0f\n', ...
    length(find(peaks.consMotif)), length(find(peaks.closeTO)), ...
    length(find(peaks.AddMotif)), length(find(ismember(peaks.motifDef, 'noMotif'))));
%% classify binding site by peak height & motif score
peak_bin = [peaks.peak_cer, peaks.peak_par];
motif_bin = [m_bin(:, 1), m_bin(:, 3)];
hasMotif_cer = motif_bin(:, 1) == 1;
hasMotif_par = motif_bin(:, 2) == 1;
peaks.hasMotif = any(motif_bin == 1, 2);
both_have_motif = all(motif_bin == 1, 2);

BSdef = cell(Npeaks, 1);
bothAboveBg = sum(peaks.isAboveBg, 2) == 2;
bothHavePeak = sum(peak_bin, 2) == 2;
if strcmp(peakStandardToUse, 'background')
    peakStandard = bothAboveBg;
elseif strcmp(peakStandardToUse, 'peak')
    peakStandard = bothHavePeak;
end

% peak in both, same motif in both
consBS_op1 = peakStandard & peaks.consMotif==1;
% peak in both, any motif in both
consBS_op2 = peakStandard & both_have_motif;
consBS = find(consBS_op1 | consBS_op2);

TOBS = find(peakStandard & peaks.closeTO==1);

addBS_cer_op1 = (peaks.AddMotif== 1 & peaks.isAboveBg(:, 1) == 1 & ...
    hasMotif_cer == 1);
addBS_cer_op2 = (peaks.closeTO == 1 & peaks.isAboveBg(:, 1) == 1 & ...
    hasMotif_cer == 1 & peaks.isAboveBg(:, 2) == 0);
addBS_cer = addBS_cer_op1 | addBS_cer_op2;

addBS_par_op1 = (peaks.AddMotif== 1 & peaks.isAboveBg(:, 2) == 1 & ...
    hasMotif_par == 1);
addBS_par_op2 = (peaks.closeTO == 1 & peaks.isAboveBg(:, 2) == 1 & ...
    hasMotif_par == 1 & peaks.isAboveBg(:, 1) == 0);
addBS_par = addBS_par_op1 | addBS_par_op2;

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

fprintf('Conserved BS: %.0f\nTurnover BS: %.0f\nAddition BS: %.0f\nNo defition: %.0f\n', ...
    length(find(ismember(BSdef, 'cons'))), ...
    length(find(ismember(BSdef, 'TO'))), ...
    length(find(ismember(BSdef, 'Addition'))), ...
    length(find(ismember(BSdef, 'noDef'))));
%% classify promoters

uniqPrs = unique(peaks.gene(peaks.hasMotif), 'stable');
n_prs = length(uniqPrs);
prT = cell(n_prs, 2);
prT(:, 1) = uniqPrs;
for i = 1:n_prs
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
pr_types = {'cons', 'TO', 'add', 'addOnly', 'noDef'};

prT.des = prT.type;

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

figure; 
spRange = [1 4];
Ng = length(pr_types);

% sum over promoter:
idxTF = find(ismember(factors, currTF));
sop = [ds_av.cer.sum_over_promoter(:, idxTF), ds_av.par.sum_over_promoter(:, idxTF)];
sop_max = max(sop, [], 2);
% sop_max = log2(sop_max+1);
gidx = ismember_smart(allGenes, prT.pr);
prT.sop_max = sop_max(gidx);
pr_groups = cell(1, length(pr_types));
for i = 1:Ng
    currGidx = find(ismember(prT.des, pr_types{i}));
    pr_groups{i} = prT.sop_max(currGidx);
end
subplot(spRange(1), spRange(2), 1);
plotSpread(pr_groups);
set(gca, 'xtick', 1:length(pr_types), 'xticklabels', pr_types);
set(gca, 'yscale', 'log');
ylabel('sum signal over promoter (max allele)');
axis square;

subplot(spRange(1), spRange(2), 2); hold on;
peakTypes = {'consMotif', 'closeTO', 'AddMotif'};
peakTypeN = sum([peaks.consMotif, peaks.closeTO, peaks.AddMotif])';
bar(peakTypeN); 
set(gca, 'xtick', 1:length(peakTypes), 'xticklabels', peakTypes);
axis square;
ylabel('# peaks');

subplot(spRange(1), spRange(2), 3); hold on;
prTypeN = nan(1, length(pr_types));
for i = 1:Ng
    prTypeN(i) = length(find(ismember(prT.des, pr_types{i})));
end
bar(prTypeN);
set(gca, 'xtick', 1:length(pr_types), 'xticklabels', pr_types);
axis square;
ylabel('# promoters');

signalBothAlleles = sum(sop, 2);
prT.sop_both = signalBothAlleles(gidx);
sum_group = nan(1, length(pr_types));
for i = 1:Ng
    currGidx = find(ismember(prT.des, pr_types{i}));
    if norm_to_total
        sum_group(i) = sum(prT.sop_both(currGidx));
    end
    if norm_to_max
        sum_group(i) = sum(prT.sop_max(currGidx));
    end
end
subplot(spRange(1), spRange(2), 4);
if norm_to_total
    prcOfTotalSignal = sum_group ./ sum(prT.sop_both)*100;
end
if norm_to_max
    prcOfTotalSignal = sum_group ./ sum(prT.sop_max)*100;
end
bar(prcOfTotalSignal);
set(gca, 'xtick', 1:length(pr_types), 'xticklabels', pr_types);
axis square;
ylabel('% of signal (linear)');

suptitle(currTF);
set(gcf, 'position', [1          41        1600         783], 'color', 'w');

tabulate(prT.TOtype(strcmp(prT.des, 'TO')))
%% peak level per group
maxPeakSum = cell(1, Ng);
groupID = {};
for i = 1:Ng
    currGroup = pr_types{i};
    idx = find(strcmp(prT.type, currGroup));
    d = vertcat(prT.maxPeakSum{idx});
    maxPeakSum{i} = d;
    groupID = vertcat(groupID, repmat({currGroup}, length(d), 1));
end
medianPeakLevel = cellfun(@nanmedian, maxPeakSum);
y = vertcat(maxPeakSum{:});
peakLevelTable = table(y, groupID, 'variableNames', {'maxSum', 'BSdef'});

%%

if flt_peaks
    savedir = [homeDir, 'checSeq_project\analyze\new_version\promoter_classification\classified_peaks_',...
        num2str(prc_top_peaks), '/'];
else
    savedir = [homeDir, 'checSeq_project\analyze\new_version\promoter_classification\classified_peaks/'];
end
save([savedir, 'peaks_', currTF, '.mat'], 'peaks');
save([savedir, 'prT_', currTF, '.mat'], 'prT');