%% allocate data structures
cd([homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\adjacent_peaks\']);
figureFolder = [pwd, '/pr_type_figs/'];
prVis_folder = [pwd, '/pr_type_prVis/'];
groups = {'cons', 'TO', 'add', 'addOnly', 'noDef'};

factors3 = {'ACE2', 'MBP1', 'REB1', 'SUT1', 'PHD1', 'SWI5', 'GCR2', 'MSN2', 'TBF1', 'SWI4', ...
    'MIG1', 'RAP1', 'PHO2', 'FKH2', 'SKN7', 'ASH1', 'PHO4'};

% NTF = length(factors3);
prTypeS = struct;
stats_TOloop = struct;
stats_TOloop.peakType = nan(NTF, 3);
stats_TOloop.prType = nan(NTF, length(groups));
stats_TOloop.prType_prc_of_signal = nan(NTF, length(groups));
peakLevels = struct;
corrOrth = struct;
for i = 1:length(groups)
    corrOrth.(groups{i}) = nan(NTF, 2);
end

%% loop over TF, run quantTO
for IDX = 1:NTF
    currTF = factors2{IDX};
    
    quant_TO_3;
    
    prTypeS.(currTF) = prT;
    set(gcf, 'invertHardCopy', 'off');
    currTitle = ['s_', num2str(IDX), '_', currTF];
    saveas(gcf, [figureFolder, currTitle, '.png'], 'png');  close(gcf);
    stats_TOloop.peakType(IDX, :) = peakTypeN;
    stats_TOloop.prType(IDX, :) = prTypeN;
    stats_TOloop.prType_prc_of_signal(IDX, :) = prcOfTotalSignal;
    stats_TOloop.medianPeakLevel(IDX, :) = medianPeakLevel;
    peakLevels.(currTF) = peakLevelTable;
    
    Npr = length(prT.pr);
    ordered_by_type = [];
    sortdClusters1 = [];
    clustering_promoters;
    for i = 1:length(groups)
        currGroup = groups{i};
        genes = prT.pr(ismember(prT.des, currGroup));
        idx = find(ismember(PrPeaks.genes, genes));
        currChangeInStrongestPeak = PrPeaks.LR(idx, basesAroundMax);
        [~, si] = sort(currChangeInStrongestPeak);
        idx_to_take = idx(si);
        ordered_by_type = [ordered_by_type; idx_to_take];
        sortdClusters1 = [sortdClusters1; repmat(i, length(idx_to_take), 1)];
    end
    si = ordered_by_type;
    plot_clustered_promoters_simple;
    currTitle = ['s_', num2str(IDX), '_', currTF];
    set(gcf, 'position', [1          41        1600         783]);
    saveas(gcf, [prVis_folder, currTitle, '.png'], 'png');  close(gcf);
    
    % correlation b/w orthologs: promoters, peaks
    for i = 1:length(groups)
        currGroup = groups{i};
        prs = prTypeS.(currTF).pr(ismember(prTypeS.(currTF).des, currGroup));
        if isempty(prs); continue; end
        otherPrs =  prTypeS.(currTF).pr(~ismember(prTypeS.(currTF).des, currGroup));
        currG = prs;
        idxPeaks = find(ismember(peaks.gene, currG) & peaks.hasMotif_sum > 0);
        peakCorr = corr(peaks.sum_cer(idxPeaks), peaks.sum_par(idxPeaks), 'rows', 'pairwise');
        prSum = nan(length(currG), 2);
        for j = 1:length(currG)
            idx = find(ismember(peaks.gene, currG{j}) & peaks.hasMotif_sum > 0);
            prSum(j, :) = [sum(peaks.sum_cer(idx)), sum(peaks.sum_par(idx))];
        end
        prCorr = corr(prSum(:, 1), prSum(:, 2), 'rows', 'pairwise');
        corrOrth.(currGroup)(IDX, :) = [prCorr, peakCorr];
    end
end
disp('done')
save('stats_TOloop.mat', 'stats_TOloop');
