% manage experiments: new 20210830

%% read bedGraph into per-base covarge maps, separate cer and par genomes
addpath(genpath(homeDir));

bedgraph_dir_path = [homeDir, 'checSeq_project\bedGraphs\'];
datadir = 'data_in\';

sp = {'cer', 'par'};

bedgraph_dir = dir(bedgraph_dir_path);
bgfiles = {bedgraph_dir.name};
bgfiles = bgfiles(contains(bgfiles, '.bedgraph'));
sampleNames = cellfun(@(x) x(1:end-9), bgfiles, 'UniformOutput', false);
Nsamples = length(bgfiles);

GenomeUsed = [datadir, 'combinedNew_wMito.fasta'];
fa = fastaread(GenomeUsed);
chrLengths = cellfun(@length, {fa.Sequence});
chrLengths_cumsum = cumsum(chrLengths);
chrLengths_cumsum_1 = [1, chrLengths_cumsum(1:end-1)+1];
load([datadir, 'varRegions_R64.mat']);
nucChr = 1:16;
chrBySp = {1:17, 18:34};

ds_raw = struct;
for i = 1:Nsamples
    disp(sampleNames{i});
    currBedGraphFile = [bedgraph_dir_path, bgfiles{i}];
    rawData = bedGraph2out(currBedGraphFile, chrLengths);
    for j = 1:size(varRegions, 1)
        currRegion = varRegions(j, :);
        rawData{currRegion(1)}(currRegion(2) : currRegion(3)) = 0;
    end
    ds_raw.cer.(sampleNames{i}) = rawData(chrBySp{1});
    ds_raw.par.(sampleNames{i}) = rawData(chrBySp{2});
end
%% normalize (only nuclear chromosomes)

multFact = 10000000;
ds_norm = struct;
Nreads = nan(Nsamples, 5);

for i = 1:Nsamples
    for j = 1:length(sp)
        currSp = sp{j};
        currRaw = ds_raw.(currSp).(sampleNames{i});
        oneVec = [currRaw{nucChr}];
        normSignal = oneVec ./ sum(oneVec)*multFact;
        signalChr = cell(1, 16);
        curr_chrLengths_cumsum = cumsum(chrLengths(chrBySp{j}));
        curr_chrLengths_cumsum_1 = [1, curr_chrLengths_cumsum(1:end-1)+1];
        for k = 1:16
            currChrLims = curr_chrLengths_cumsum_1(k):curr_chrLengths_cumsum(k);
            signalChr{k} = normSignal(currChrLims);
        end
        ds_norm.(currSp).(sampleNames{i}) = signalChr;
    end
    curr_totalReads = sum([sum(cellfun(@sum, ds_raw.cer.(sampleNames{i}))), sum(cellfun(@sum, ds_raw.par.(sampleNames{i})))]);
    curr_mito = sum([sum(ds_raw.cer.(sampleNames{i}){17}), sum(ds_raw.par.(sampleNames{i}){17})]);
    curr_notmito = curr_totalReads - curr_mito;
    curr_cer_not_mito = sum(cellfun(@sum, ds_raw.cer.(sampleNames{i})(1:16)));
    curr_par_not_mito = sum(cellfun(@sum, ds_raw.par.(sampleNames{i})(1:16)));
    Nreads(i, :) = [curr_totalReads, curr_mito, curr_notmito, curr_cer_not_mito, curr_par_not_mito];
end
        
NreadsTable = array2table(Nreads, 'variableNames', {'total', 'mito', 'not_mito', 'cer_not_mito', 'par_not_mito'}, ...
    'RowNames', sampleNames);
save('ds_norm.mat', 'ds_norm', '-v7.3');
%% average repeats

bioSampleT = readtable([datadir, 'BioSample.xlsx']);
idx = find(ismember(bioSampleT.sample_name, sampleNames));
bioSampleT = bioSampleT(idx, :);
factors = unique(bioSampleT.transcription_factor, 'stable');
NTF = length(factors);
strainList = struct;
for i = 1:NTF
    strainList.(factors{i}) = bioSampleT.sample_name(strcmp(bioSampleT.transcription_factor, factors{i}))';
end

av_norm = struct;
genomeLen = [sum(chrLengths(1:16)); ...
   sum(chrLengths(18:33))];
chrs = {1:16, 18:33};
for i = 1:length(factors)
    currSamples = strainList.(factors{i});
    for j = 1:2
        currVec = zeros(length(currSamples), genomeLen(j));
        for k = 1:length(currSamples)
            currVec(k, :) = [ds_norm.(sp{j}).(currSamples{k}){:}];
        end
        meanVec = nanmean(currVec, 1);
        
        currChr = chrs{j};
        curr_chrLengths_cumsum = cumsum(chrLengths(currChr));
        curr_chrLengths_cumsum_1 = [1, curr_chrLengths_cumsum(1:end-1)+1];
        signalChr = cell(1, 16);
        for k = 1:16
            currChrLims = curr_chrLengths_cumsum_1(k):curr_chrLengths_cumsum(k);
            signalChr{k} = meanVec(currChrLims);
        end
        
        av_norm.(sp{j}).norm.(factors{i}) = signalChr;
    end
end
    
save('av_norm.mat', 'av_norm', '-v7.3');
%% sum over promoter

prLength = 700;

global tss_struct;
load([datadir, 'tss_struct.mat']);
ds_av = struct;
for i = 1:length(sp)
    currSp = sp{i};
    curr_sop = sumOnPro(av_norm.(currSp), prLength, currSp);
    curr_sop = cell2mat(struct2cell(curr_sop))';
    ds_av.(currSp).sum_over_promoter = curr_sop;
end
   
save('ds_av.mat', 'ds_av');
%% organize full hybrid genome

fa2 = fa([1:16, 18:33]);
% genome
genomes.hyb = {fa2.Sequence}';

% TSS
tss_hyb = nan(6701*2, 3);
tss_hyb(1:6701, :) = tss_struct.cer;
tss_hyb(6702:13402, :) = tss_struct.par;
tss_hyb(6702:13402, 1) = tss_hyb(6702:13402, 1) + 16;

% intergenic region
inter_hyb = [tss_struct.intergenic.cer; tss_struct.intergenic.par];

tss_struct.hyb = tss_hyb;
tss_struct.intergenic.hyb = inter_hyb;

%% 7-mer motif score
motifStruct = struct;
nmer = 7;
currMer = ['mer', num2str(nmer)];
factors = fields(av_norm.cer.norm);
profile = nan(sum(chrLengths([1:16, 18:33])), length(factors));
for i = 1:length(factors)
    profile(:, i) = [av_norm.cer.norm.(factors{i}){:}, av_norm.par.norm.(factors{i}){:}]';
end
out = mer_occupancy_general_wholeHyb(profile, nmer);
% clear profile;
for i = 1:length(factors)
    motifStruct.hyb.(currMer).meanSignalPerMotif.(factors{i}) = out.score(:, i);
end
motifStruct.hyb.(currMer).motifs_seq = out.mers.seq;
motifStruct.hyb.(currMer).motifs_seqRC = out.mers.rcSeq;
motifStruct.hyb.(currMer).table = out.mers;
motifStruct.hyb.(currMer).genomeXmers = out.sc4nmer;
motifStruct.hyb.promotersRegion = out.promotersRegion;
%
%% save PWMs

seqlogo_dir = 'seqlogos\';
if ~exist(seqlogo_dir, 'dir')
    mkdir(seqlogo_dir);
end
motNum = 20; % top 7-mer motif sequences to include in the seqlogo
for i = 1:length(factors)
    currTF = factors{i};
    disp(currTF);
    [finalFreq,totAGCT,Positions2Plot] = seqlogo_wrap(motifStruct.hyb, currTF, motNum, currMer);
    motifStruct.hyb.(currMer).PFM.(currTF) = finalFreq(:, Positions2Plot);
    set(gcf, 'invertHardCopy', 'off');
    saveas(gcf, [seqlogo_dir, currTF, '.png'], 'png');  close(gcf);
end
save('motifStruct.mat', 'motifStruct');  
%% Align orthologous promoters

S288c_R64 = fastaread([datadir, 'S288C_reference_sequence_R64-1-1_20110203.fasta']);
CBS432 = fastaread([datadir, 'CBS432.genome.fa']);

genomes.S288c_R64_cer = {S288c_R64.Sequence};
genomes.S288c_R64_cer = genomes.S288c_R64_cer(1:16)';
genomes.CBS432_par = {CBS432.Sequence};
genomes.CBS432_par = genomes.CBS432_par(1:16)';

load([datadir, 'GP.mat']);

signal_on_aligned_promoters_part1;
signal_on_aligned_promoters_part2;

%% Signal around motif (potential binding sites, PBS)

pbs3_1_sortSignal_ver2;
pbs3_2_seqConservationAtCore_ver2;
pbs3_3_signalAroundMotif_ver2;

%% generate aligned motif score data structure: alMotifScoresS

score_motif_by_kmers_loop;

%% Generate peak table
currTF = 'REB1'; % for example
MinPeakHeightFixed = 594;
peaks_vs_motifScore;

%% mutations in binding sites

effect_of_fixated_mut_loop;

%% promoter classification

quant_TO_loop_plot;