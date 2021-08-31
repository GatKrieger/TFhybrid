%% signal on aligned promoters

averageRepeats = 1;

if ~exist('alignedLoc'); load('alignedLoc_1_3.mat'); end
if ~exist('alPromSeq'); load('alPromSeq_1_3.mat'); end
global alignedLocV alignedLocGene signalAligned
alignedLocV = struct;

%% take once the promoter of head2head genes
alignedLoc2 = alignedLoc;
genesWTss = find(~isnan(tss_struct.cer(:, 1)) & ~isnan(tss_struct.par(:, 1)));
if ~exist('head2head'); run('type_of_promoter_tail2tail_etc.m'); end
head2head_WTss = intersect(genesWTss, head2head);
length(head2head_WTss);
for i = 1:length(head2head)
    g1 = find(genesWTss == head2head(i));
    g2 = find(genesWTss == head2head_neighbor(i));
    if ~isempty(g1) & ~isempty(g2)
        alignedLoc2.cer{g1} = [];
        alignedLoc2.par{g1} = [];
        alPromSeq(g1).Alignment = [];
    end
end
c = [alignedLoc2.cer{:}];
p = [alignedLoc2.par{:}];
alignedLocV.cer = c(1, :);
alignedLocV.par = p(1, :);
alignedLocGene = c(2, :);
%% unique indices
wholeAlSeq = [alPromSeq.Alignment];
idxNotNan = find(~isnan(alignedLocV.cer));
length(idxNotNan);
[uniqueC, ia, ic] = unique(alignedLocV.cer(idxNotNan), 'stable');
uniqueP = alignedLocV.par(idxNotNan);
uniqueP = uniqueP(ia);

vcer = nan(1, length(wholeAlSeq));
vcer(idxNotNan(ia)) = uniqueC;
vpar = nan(1, length(wholeAlSeq));
vpar(idxNotNan(ia)) = uniqueP;

idxNanVpar = find(isnan(vpar));
vpar(idxNanVpar) = 0;
% correction:
for i = 1:length(vpar)
    if vpar(i) == 0
        if vpar(i-1) == 0; continue; end
        n = 1;
        while vpar(i + n) == 0 & n+i < length(vpar)
            n = n+1;
        end
        if abs(vpar(i + n) - vpar(i-1)) > 1
            if vpar(i + n) - vpar(i-1) > 1
                if length(vpar(i-1:i+n)) == length(vpar(i-1) : vpar(i+n))
                    vpar(i-1:i+n) = vpar(i-1) : vpar(i+n);
                end
            elseif vpar(i + n) - vpar(i-1) < -1
                if length(vpar(i-1:i+n)) == length(vpar(i-1) :-1: vpar(i+n))
                    vpar(i-1:i+n) = vpar(i-1) :-1: vpar(i+n);
                end
            end
        end
    end
end
vpar(vpar == 0) = nan;

alignedLocV.cer = vcer;
alignedLocV.par = vpar;
clear vcer vpar;
fprintf('size of aligned promoters: %.0f\n', length(alignedLocV.cer));
% alg = nan(1, length(wholeAlSeq));
% alg(idxNotNan(ia)) = alignedLocGene(idxNotNan(ia));
%% align ChEC signal
if ~exist('ds_norm', 'var'); load('ds_norm.mat'); end
smoothingSize = 20;
slidingWindow = ones(smoothingSize, 1);
l = length(alignedLocGene);

alignedPrLen = length(alignedLocGene);
genLen = [sum(chrLength.S288c_R64_cer), sum(chrLength.CBS432_par)];

if averageRepeats
    wantedSamples = factors;
else
    wantedSamples = struct2array(strainList);
end
% F = fields(ds_norm.cer);
signalAligned = struct;
signalAlignedSm = struct;
for i = 1:length(sp)
    signalAligned.(sp{i}) = nan(l, length(wantedSamples));
    signalAlignedSm.(sp{i}) = nan(l, length(wantedSamples));
    v = alignedLocV.(sp{i});
    idxNotNan = find(~isnan(v));
    for j = 1:length(wantedSamples)
        if averageRepeats
            idxReps = strainList.(wantedSamples{j});
%             idxReps = F(contains(F, wantedSamples{j}));
            origSignal = nan(length(idxReps), genLen(i));
            for k = 1:length(idxReps)
                origSignal(k, :) = [ds_norm.(sp{i}).(idxReps{k}){:}];
            end
            origSignal = nanmean(origSignal, 1);
        else
            origSignal = [ds_norm.(sp{i}).(wantedSamples{j}){:}];
        end
        currSig = nan(1, l);
        currSig(idxNotNan) = origSignal(v(idxNotNan));
        signalAligned.(sp{i})(:, j) = currSig;
        
        % smooth data with zeros, not nans. smoothing nans doesn't work...
        currSig(isnan(currSig)) = 0;
%         smSignal = smoothdata(currSig, 'movmean', smoothingSize);
        smSignal = conv(currSig, slidingWindow, 'same');
        signalAlignedSm.(sp{i})(:, j) = smSignal;
    end
end
%%
save('alignedLocV.mat', 'alignedLocV', '-v7.3');
save('alignedLocGene.mat', 'alignedLocGene');
save('alPromSeq_1_3.mat', 'alPromSeq');
if averageRepeats == 0
    save('signalAligned_reps.mat', 'signalAligned', '-v7.3');
    save('signalAlignedSm_reps.mat', 'signalAlignedSm', '-v7.3');
    signalAligned_reps_samples = wantedSamples;
    save('signalAligned_reps_samples.mat', 'signalAligned_reps_samples');
else
    save('signalAligned.mat', 'signalAligned', '-v7.3');
    save('signalAlignedSm.mat', 'signalAlignedSm', '-v7.3');
end
%% aligned sequence

wholeAlignment = horzcat(alPromSeq(:).Alignment);
alignedSeq = struct;
alignedSeq.cer = wholeAlignment(1, :);
alignedSeq.par = wholeAlignment(3, :);
save('alignedSeq.mat', 'alignedSeq');
%% nucleosomes on aligned Promoters
if ~exist('nucRealign'); load([homeDir, 'NucleosomesRaw/nucRealign.mat']); end
global nucAligned
nucAligned = struct;
longName_alsoHyb = {'S288c_R64_cer', 'CBS432_par', 'S288c_R64_cer', 'CBS432_par'};
sp_alsoHyb = {'cer', 'par', 'cer', 'par'};
divisionFac = 80;
for i = 1:length(gb)
    currGenome = genomes.(longName_alsoHyb{i});
    chrLen = cellfun(@length, currGenome);
    v = alignedLocV.(sp_alsoHyb{i});
    idxNotNan = find(~isnan(v));
    origSignal = [nucRealign.(['profileNorm_', gb{i}]){:}];
    currSig = nan(1, l);
    currSig(idxNotNan) = origSignal(v(idxNotNan));
    currSig = currSig./divisionFac;
    nucAligned.(gb{i}) = currSig;
end
load([homeDir, 'Tirosh data/nucFridman_allgenome.mat']);
currGenome = genomes.S288c_R64_cer;
chrLen = cellfun(@length, currGenome);
v = alignedLocV.cer;
idxNotNan = find(~isnan(v));
origSignal = [inputdataMedian{:}];
currSig = nan(1, l);
currSig(idxNotNan) = origSignal(v(idxNotNan));
nucAligned.cer2 = currSig./divisionFac;
save('nucAligned.mat', 'nucAligned');
clear inputdataMedian;
%% UCSC
ucsccons = load([homeDir, 'Ext_data/UCSC_genome_browser_conservation_track/data.mat']);
currGenome = genomes.S288c_R64_cer;
chrLen = cellfun(@length, currGenome);
v = alignedLocV.cer;
idxNotNan = find(~isnan(v));
origSignal = [];
for i = 1:16
    origSignal = [origSignal, ucsccons.profile{i}'];
end
currSig = nan(1, l);
currSig(idxNotNan) = origSignal(v(idxNotNan));
global ucsccons_aligned;
ucsccons_aligned = currSig;
save('ucsccons_aligned.mat', 'ucsccons_aligned');
clear v origSignal idxNotNan