% plot mutation frequency and changes in binding to all factors

wantedMinScoreTh = -1;
addition = [];

dataDir = [homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\binding_site_mutations\'...
    'mat_files_onlyOneMut_inVitro_std', stdStr, '/'];
load([dataDir, 'absSumPerPosMutAll-1.mat']);
load([dataDir, 'LRperPosMutAll-1.mat']);

flankingBases = (length(fields(LRperPosMutAll.REB1)) - 7)/2;
factors2 = factors;
factors2(ismember(factors2, {'RPB9', 'FHL1', 'STB3', 'PHO2', 'HAP4', 'HMS2', 'FKH2', 'TEC1', 'MIG2'})) = [];
NTF = length(factors2);

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

Nalt_per_nuc = struct;
deltaBind = struct;
mutTypes = {'nonFixed_same', 'nonFixed_bw', 'fixed_bw'};
for i = 1:length(mutTypes)
    Nalt_per_nuc.(mutTypes{i}) =  cell(NTF, 1);
    deltaBind.(mutTypes{i}) =  cell(NTF, 1);
end
maxLen2 = maxLen + flankingBases*2;

for i = 1:NTF
    currTF = factors2{i};
    currF = LRperPosMutAll.(currTF);
    f = fields(currF);
    Nf = length(f);
    deltaBind.nonFixed_same{i} = nan(4, maxLen2);
    deltaBind.nonFixed_bw{i} = nan(4, maxLen2);
    deltaBind.fixed_bw{i} = nan(4, maxLen2);
    
    Nalt_per_nuc.nonFixed_same{i} = nan(4, maxLen2);
    Nalt_per_nuc.nonFixed_bw{i} = nan(4, maxLen2);
    Nalt_per_nuc.fixed_bw{i} = nan(4, maxLen2);
    
    currMotif = PFMtoUse.(currTF);
    gapSize = maxLen2 - Nf;
    halfGapSize = floor(gapSize / 2);
    posToFill = halfGapSize+1 : halfGapSize+Nf;
    
    for j = 1:length(f)
        cPos = f{j};
        currPosToFill = posToFill(j);
        % non-fixed mutation, same site
        currPos = LRperPosMutAll.(currTF).(cPos).bwSp;
        NsitesPerNuc = cellfun(@length, currPos);
        NsitesPerNuc(1) = []; % gap
        currCons = consensusSeq.consNum.(currTF)(j);
        if isnan(currCons)
            [~, currCons] = max(NsitesPerNuc);
        end
        %         NsitesPerNuc(currCons) = nan;
        Nalt_per_nuc.nonFixed_same{i}(:, currPosToFill) = NsitesPerNuc;
        currPos = LRperPosMutAll.(currTF).(cPos).bwSp(2:5);
        deltaBind.nonFixed_same{i}(:, currPosToFill) = cellfun(@mean, currPos);
        
        % fixed mutation, between sites
        currPos = absSumPerPosMutAll.(currTF).(cPos).inGen(2:5);
        NsitesPerNuc = cellfun(@length, currPos);
        if isnan(currCons)
            [~, currCons] = max(NsitesPerNuc);
        end
        %         NsitesPerNuc(currCons) = nan;
        Nalt_per_nuc.fixed_bw{i}(:, currPosToFill) = NsitesPerNuc;
        deltaBind.fixed_bw{i}(:, currPosToFill) = ...
            [absSumPerPosMutAll.(currTF).(cPos).inGenConsLR{2:5}];
        
        % non-fixed mutation, between sites
        currPos = absSumPerPosMutAll.(currTF).(cPos).bwSp(2:5);
        NsitesPerNuc = cellfun(@length, currPos);
        if isnan(currCons)
            [~, currCons] = max(NsitesPerNuc);
        end
        %         NsitesPerNuc(currCons) = nan;
        Nalt_per_nuc.nonFixed_bw{i}(:, currPosToFill) = NsitesPerNuc;
        mean4 = cellfun(@mean, currPos);
        delta4 = mean4 - mean4(currCons);
        deltaBind.nonFixed_bw{i}(:, currPosToFill) = delta4;
    end
end
%% plot

figure;
spRange = [1 6];
n=1;
ax = cell(1,6);
cmDelta = flipud(cbrewer2('BrBg'));
bgColor = [.7 .7 .7];
deltas = -10:5:10;
xticks = ceil(maxLen/2) + deltas;
xlabels = cellfun(@num2str, num2cell(deltas), 'uniformoutput', false);

% red squares at the consensus sequence
consMat = struct;
for i = 1:NTF
    currTF = factors2{i};
    currCons = consensusSeq.consNum.(currTF);
    currLen = length(currCons);
    idxes = sub2ind([4 currLen], currCons, 1:currLen);
    idxes(isnan(idxes)) = [];
    currMat = zeros(4, currLen);
    currMat(idxes) = 1;
    gapSize = maxLen - currLen;
    halfGapSize = floor(gapSize / 2);
    posToFill = halfGapSize+1 : halfGapSize+currLen;
    currMat2 = zeros(4, maxLen);
    currMat2(:, posToFill) = currMat;
    consMats.(currTF) = currMat2;
end
seqLengths = structfun(@length, consensusSeq.consNum);
maxLen = max(seqLengths);
consMat = cell2mat(struct2cell(consMats));
idxes = find(consMat);
[row, col] = ind2sub(size(consMat), idxes);
idxes2 = [row, col];

ttls = {'non-fixed, same', 'non-fixed, between', 'fixed, between'};
for i = 1:length(mutTypes)
    s = subplot(spRange(1), spRange(2), n); n=n+1;
    m = cell2mat(deltaBind.(mutTypes{i}));
    Ns = cell2mat(Nalt_per_nuc.(mutTypes{i}));
    ii = imagesc(m); hold on;
    set(ii, 'alphadata', Ns >= minSites); set(gca, 'color', bgColor);
%     cb = colorbar('southoutside');
%     cb.Label.String = 'log2 (alt / cons)'; cb.FontSize = fs;
    title({ttls{i}, '\Deltabinding'});
    colormap(s, cmDelta);
    caxis([-5 5]);
    for iii = 1:NTF; rectangle('Position', rectPos4f(iii, :)); end
    ytickPos = 2.5:Ncats:NTF*Ncats; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
    set(gca, 'xtick', xticks, 'xticklabels', xlabels);
    
    for iii = 1:length(idxes2)
        rectangle('Position', [idxes2(iii, 2)-.5, idxes2(iii, 1)-.5, 1, 1], 'edgeColor', 'r');
    end
    
    s = subplot(spRange(1), spRange(2), n); n=n+1;
    m = Ns;
    ii = imagesc(m); hold on;
    set(ii, 'alphadata', Ns >= minSites); set(gca, 'color', bgColor);
%     cb = colorbar('southoutside');
%     cb.Label.String = '# sites'; cb.FontSize = fs;
    title({ttls{i}, '# sites'});
    caxis([3 2000])
    set(gca, 'ColorScale', 'log');
    for iii = 1:NTF; rectangle('Position', rectPos4f(iii, :)); end
    ytickPos = 2.5:Ncats:NTF*Ncats; set(gca, 'ytick', ytickPos, 'yticklabels', factors2);
    set(gca, 'xtick', xticks, 'xticklabels', xlabels);
end