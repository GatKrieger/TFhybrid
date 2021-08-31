%% log-ratio of alt/cons: in fixed and in non-fixed
Ncats = 2;
% reorderedfactors2 = {'REB1', 'MBP1', 'TBF1','MSN2','RAP1', 'MIG1', 'SUT1', 'ACE2', 'SWI4', 'SKN7', 'MIG2', ...
%     'SWI5', 'ABF1', 'SOK2', 'TOD6', 'PHD1', 'FKH1', 'PHO4', 'COM2', 'ASH1', 'SUT2', 'GCR1', 'RGT1', 'DOT6', ...
%     'GCR2', 'YAP2', 'TEC1', 'SFP1'};
reorderedfactors2 = {'ABF1','REB1','TOD6','RAP1','TEC1','FKH1','SUT1','DOT6','SWI5','MIG1','MIG2','MSN2','SWI4','GCR1',...
    'SKN7','SOK2','FKH2','RGT1','ACE2','PHO4','GCR2','STB3','PHD1','TBF1','HAP4','MBP1','FHL1','ASH1','PHO2'};

dataDir = [homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\binding_site_mutations\'...
    'mat_files_onlyOneMut_inVitro_std1/'];
load([dataDir, 'absSumPerPosMutAll-1.mat']);
load([dataDir, 'LRperPosMutAll-1.mat']);
load([dataDir, 'rawDataAll-1.mat']);

NTF = length(reorderedfactors2);
seqLengths = structfun(@length, PFMtoUse);
maxLen = max(seqLengths) + flankingBases*2;
cvals = cell(NTF, 1);
cP = cell(NTF, 1);
rectPos = nan(NTF, 4);
NsitesPerPos = cell(NTF, 1);
whatComparison = 'bwSites';
% whatComparison = 'wiSite';


for i = 1:NTF
    currTF = reorderedfactors2{i};
    if ~isfield(absSumPerPosMutAll, currTF)
        cvals{i} = nan(Ncats, maxLen);
        cP{i} = nan(Ncats, maxLen);
        NsitesPerPos{i} = nan(1, maxLen);
        continue;
    end
    Npos = length(fields(absSumPerPosMutAll.(currTF)));
    deltaMean = nan(2, Npos);
    deltaPvals = nan(2, Npos);
    NsitesInPosNonfixed = nan(Npos, 1);
    NsitesInPosFixed = nan(Npos, 1);
    for j = 1:Npos
        currPos = ['pos', num2str(j)];
        deltaMean(2, j) = absSumPerPosMutAll.(currTF).(currPos).inGenConsLR_allAlt;
        deltaPvals(2, j) = absSumPerPosMutAll.(currTF).(currPos).inGenConsPval_allAlt;
        NsitesInPosFixed(j) = sum(cellfun(@length, absSumPerPosMutAll.(currTF).(currPos).inGen))/2;
        switch whatComparison
            case 'bwSites'
                deltaMean(1, j) = absSumPerPosMutAll.(currTF).(currPos).bwSpConsLR_allAlt;
                deltaPvals(1, j) = absSumPerPosMutAll.(currTF).(currPos).bwSpConsPval_allAlt;
            case 'wiSite'
                deltaMean(1, j) = LRperPosMutAll.(currTF).(currPos).bwSpConsLR_allAlt;
                deltaPvals(1, j) = LRperPosMutAll.(currTF).(currPos).bwSpConsPval_allAlt;
        end
        NsitesInPosNonfixed(j) = sum(cellfun(@length, LRperPosMutAll.(currTF).(currPos).bwSp));
    end
    m2 = nan(Ncats, maxLen);
    sizeDiff = maxLen - Npos;
    toGap = round(sizeDiff/2);
    m2(:, toGap+1: toGap+ Npos) = deltaMean;
    cvals{i} = m2;
    currMotifLength = length(PFMtoUse.(currTF));
    rectPos(i, :) = [toGap+flankingBases+.5, .5+i*Ncats - Ncats, currMotifLength, Ncats];
    
    m2 = nan(Ncats, maxLen);  m2(:, toGap+1: toGap+ Npos) = deltaPvals; cP{i} = m2;
    m2 = nan(2, maxLen);  
    m2(1, toGap+1: toGap+ Npos) = NsitesInPosNonfixed';
    m2(2, toGap+1: toGap+ Npos) = NsitesInPosFixed';
    NsitesPerPos{i} = m2;
end

%%
deltas = -10:5:10;
xticks = ceil(maxLen/2) + deltas;
xlabels = cellfun(@num2str, num2cell(deltas), 'uniformoutput', false);

m = cell2mat(cvals);
figure;
tl = tiledlayout(1, 3);
ax1 = nexttile();
i1 = imagesc(m); hold on;
set(i1, 'alphadata', ~isnan(m));
cm2 = cbrewer2('RdBu'); cm2 = cm2(end:-1:1, :);
colormap(ax1, cm); set(gca, 'color', [.7 .7 .7]);
set(gca, 'ytick', 1.5:Ncats:NTF*Ncats, 'yticklabels', reorderedfactors2);
caxis([-5 5]);
plot(xlim, [2.5:Ncats:NTF*Ncats; 2.5:Ncats:NTF*Ncats]', 'k');
plot(xlim, [1.5:Ncats:NTF*Ncats; 1.5:Ncats:NTF*Ncats]', '-', 'color', [.7 .7 .7]);
for i = 1:NTF
    if ~all(isnan(rectPos(i, :)))
        rectangle('Position', rectPos(i, :));
    end
end
cb = colorbar('southoutside'); cb.Label.String = 'alternative / consensus log2 ratio';
set(gca, 'xtick', xticks, 'xticklabels', xlabels);
% xlabel('Distance from motif center (bp)');

mPval = cell2mat(cP);
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(mPval);

ax2 = nexttile();
m = -log10(adj_p);
i2 = imagesc(m); hold on;
set(i2, 'alphadata', ~isnan(m));
set(gca, 'color', [.7 .7 .7]);
cb = colorbar('southoutside'); cb.Label.String = '-log10(p-value)';
set(gca, 'ytick', 1.5:Ncats:NTF*Ncats, 'yticklabels', []);
plot(xlim, [2.5:Ncats:NTF*Ncats; 2.5:Ncats:NTF*Ncats]', 'k');
plot(xlim, [1.5:Ncats:NTF*Ncats; 1.5:Ncats:NTF*Ncats]', '-', 'color', [.7 .7 .7]);
set(gca, 'xtick', xticks, 'xticklabels', xlabels);
xlabel('Distance from motif center (bp)');

ax3 = nexttile();
m = cell2mat(NsitesPerPos);
m = log10(m);
m(isinf(m)) = 0;
i3 = imagesc(m); hold on;
set(i3, 'alphadata', ~isnan(m));
set(gca, 'color', [.7 .7 .7]);
cb = colorbar('southoutside'); cb.Label.String = '# log10 sites';
set(gca, 'ytick', 1.5:Ncats:NTF*Ncats, 'yticklabels', []);
plot(xlim, [2.5:Ncats:NTF*Ncats; 2.5:Ncats:NTF*Ncats]', 'k');
plot(xlim, [1.5:Ncats:NTF*Ncats; 1.5:Ncats:NTF*Ncats]', '-', 'color', [.7 .7 .7]);
set(gca, 'xtick', xticks, 'xticklabels', xlabels);
% xlabel('Distance from motif center (bp)');

set(gcf, 'color', 'w');
% linkaxes([ax1,ax2,ax3], 'y')