%% 
run_classifier = 0;
if run_classifier; quant_TO_loop1; end
%% plot number of promoters per promoter type
% and percent of binding signal per promoter type
load stats_TOloop.mat
prcSignal = stats_TOloop.prType_prc_of_signal;
[~, si] = sort(prcSignal(:, 1));
colors = lines(length(groups));
colors = [96, 166, 213;...
    138, 212, 244;...
    173, 123, 183;...
    232, 149, 112;...
    184, 76, 97];
colors = colors./255;
figure; 
biggerPanel = 3;
tl = tiledlayout(1, biggerPanel*2+1);

groups_long = {'Conserved', 'Turnover', 'Unbalanced', 'Fully unbalanced', 'not defined'};
ax2 = nexttile([1 biggerPanel]);
prTypePrc = stats_TOloop.prType ./ sum(stats_TOloop.prType, 2)*100;
b = barh(prTypePrc(si, :), 'stacked');
for i = 1:length(groups)
    b(i).FaceColor  = colors(i, :);
end
set(gca, 'ytick', 1:NTF, 'yticklabels', factors2(si));
% legend(groups_long, 'location', 'eastoutside')
xlabel('% promoters');
axis([0 100 .5 NTF+.5])

ax3 = nexttile([1 biggerPanel]);
b2 = barh(prcSignal(si, :), 'stacked');
set(gca, 'ytick', 1:NTF, 'yticklabels', []);
xlabel('% signal');
axis([0 100 .5 NTF+.5])
for i = 1:length(groups)
    b2(i).FaceColor  = colors(i, :);
end
legend(groups_long, 'location', 'southoutside', 'box', 'off')

ax4 = nexttile();
Npromoters = sum(stats_TOloop.prType, 2);
barh(Npromoters(si, :));
set(gca, 'xscale', 'log');
set(gca, 'ytick', 1:NTF, 'yticklabels', []);
xlim([100, 900])
ylim([1, NTF])
xlabel('# promoters')

linkaxes([ax2, ax3, ax4], 'y')

set(gcf, 'position', [ 520   157   451   641], 'color', 'w');
%% 
% total promoter signal per class
Ngrp = length(groups);
prGroupS = cell(NTF, Ngrp);
for i = 1:NTF
    currTF = factors2{i};
    for j = 1:Ngrp
        currGidx = find(ismember(prTypeS.(currTF).des, groups{j}));
        d = prTypeS.(currTF).sop_max(currGidx);
%         d = log2(d+1);
        prGroupS{i, j} = d;
    end
end
prGroupT = array2table(prGroupS, 'variableNames', groups, 'RowNames', factors2);

totalSum = cellfun(@sum, prGroupS);
SperTF = sum(totalSum, 2);
medianGroup = cellfun(@median, prGroupS);
medianGroupT = array2table(medianGroup, 'variableNames', strrep(groups_long, ' ', '_'));
writetable(medianGroupT, 'signal_per_group_raw.csv');
% figure; plotSpread(medianGroup)
% set(gca, 'yscale', 'log')

normA = cell(NTF, Ngrp);
for i = 1:NTF
    for j = 1:Ngrp
        d = prGroupS{i,j}/SperTF(i)*1E6;
        d = log2(d+1);
        normA{i,j} = d;
    end
end
medianGroup2 = cellfun(@median, normA);
medianGroup2T = array2table(medianGroup2, 'variableNames', strrep(groups_long, ' ', '_'));
writetable(medianGroup2T, 'signal_per_group.csv');

figure; 
spRange = [1 3];
subplot(spRange(1), spRange(2), 1);
violinplot(medianGroup2);
% set(gca, 'yscale', 'log');
set(gca, 'xtick', 1:length(groups), 'xticklabels', groups_long); xtickangle(45);
ylabel('log2 sum of signal on promoter');
% axis square;
% run multiple_comp.R
cld = {'b', 'c', 'c', 'ab', 'a'};
err = 0.1;
for i = 1:length(cld)
    text(i-err, 6.5, cld{i});
end
ylim([6 14]);

%% Peak signal per class
peakLevelsNormMed = nan(NTF, Ng);
for i = 1:NTF
    currTF = factors2{i};
    peakLevels.(currTF).maxSum_norm = peakLevels.(currTF).maxSum ./ sum(peakLevels.(currTF).maxSum)*1E6;
    for j = 1:Ng
        currGroup = groups{j};
        idx = find(strcmp(peakLevels.(currTF).BSdef, currGroup));
        d = peakLevels.(currTF).maxSum_norm(idx);
        d = log2(d+1);
        peakLevelsNormMed(i, j) = nanmedian(d);
    end
end

subplot(spRange(1), spRange(2), 2);
violinplot(peakLevelsNormMed(:, 1:4));
set(gca, 'xtick', 1:length(groups), 'xticklabels', groups_long(1:4)); xtickangle(45);
ylabel('log2 signal on binding site');
% [~, ~, statsAnova] = anova1(peakLevelsNormMed);
% [c,~,~,gnames] = multcompare(statsAnova);
cld = {'b', 'b', 'b', 'a'};
err = 0.1;
for i = 1:length(cld)
    text(i-err, 6.5, cld{i});
end
ylim([6 14]);

%% Number of binding sites per class
Nbs = cell(NTF, Ng);
for i = 1:NTF
    currTF = factors2{i};
    for j = 1:Ng-1
        currGroup = groups{j};
        idx = find(strcmp(prTypeS.(currTF).type, currGroup));
%         Nbs{i,j} = cellfun(@length, prTypeS.(currTF).BSdef(idx));
        Nbs{i,j} = prTypeS.(currTF).NBS(idx);
    end
end

Nbs_mean = cellfun(@mean, Nbs);
subplot(spRange(1), spRange(2), 3);
violinplot(Nbs_mean(:, 1:4));
set(gca, 'xtick', 1:length(groups), 'xticklabels', groups_long(1:4)); xtickangle(45);
ylabel('# binding sites');
% [~, ~, statsAnova] = anova1(Nbs_mean);
% [c,~,~,gnames] = multcompare(statsAnova);
cld = {'b', 'c', 'd', 'a'};
err = 0.1;
for i = 1:length(cld)
    text(i-err, 1.1, cld{i});
end
set(gcf, 'position', [520   516   732   282], 'color', 'w');

%%

allN = cell(1, Ng-1);
cats = cell(1, Ng-1);
for i = 1:Ng-1
    allN{i} = vertcat(Nbs{:, i});
    cats{i} = repmat(groups_long(i), length(allN{i}), 1);
end
allNvec = vertcat(allN{:});
catVec = vertcat(cats{:});
figure; violinplot(allNvec, catVec);

