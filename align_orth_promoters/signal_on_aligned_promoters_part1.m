% align promoters:configs
wantedRegion = 'start';
minPromLen = 400;
maxPromLen = 1500;
downStreamBases = 50;
genesWTss = find(~isnan(tss_struct.cer(:, 1)) & ~isnan(tss_struct.par(:, 1)));
spLongName = {'S288c_R64_cer', 'CBS432_par'};
alignedLoc = struct;
alignedLoc.cer = cell(length(genesWTss), 1);
alignedLoc.par = cell(length(genesWTss), 1);
wholeGen.cer = [genomes.S288c_R64_cer{:}];
wholeGen.par = [genomes.CBS432_par{:}];
alPromSeq = struct;
%%
if ~exist('genomeToChrS', 'var')
    genomeToChrS_generate
end
%% define promoter region

promDef = struct;
for j = 1:length(sp)
    currIntR = tss_struct.intergenic.(sp{j});
    currIntR(isnan(currIntR)) = minPromLen;
    currIntR(currIntR <= 0) = minPromLen;
    currIntR(currIntR > maxPromLen) = maxPromLen;
    promDef.(sp{j}) = currIntR + downStreamBases;
    
    % take 1000 bp upstream
%     promDef.(sp{j}) = maxPromLen * ones(6701, 1);
end

%%
for i = 1:length(genesWTss)
    %% pairwise sequence alignment for each promoter using matlab 'localalign'
    currGn = genesWTss(i);
    disp([num2str(currGn), ', ', allGenes{currGn}]);
    seq = cell(1,2);
    region = nan(2,3);
    regionInGenome = nan(2,2);
    strand = ones(1,2);
    for j = 1:2
        upStreamBases = promDef.(sp{j})(currGn);
        [seq{j}, loc] = extractGenomicRegion2(currGn, wantedRegion, upStreamBases, downStreamBases,...
            tss_struct.(sp{j}), genomes.(spLongName{j}), 0);
        region(j, :) = loc;
        regionInGenome(j, 1) = find(genomeToChrS.(sp{j})(2, :) == loc(1) & genomeToChrS.(sp{j})(1, :) == loc(2));
        regionInGenome(j, 2) = find(genomeToChrS.(sp{j})(2, :) == loc(1) & genomeToChrS.(sp{j})(1, :) == loc(3));
        currLoc = tss_struct.(sp{j})(currGn, :);
        if currLoc(2) > currLoc(3)
            strand(j) = 0;
        end
    end
    [score, alignment, startat] = swalign(seq{1}, seq{2}, 'gapOpen', 10, 'extendGap', .5, 'alphabet', 'NT');
    alPromSeq(i).geneID = currGn;
    alPromSeq(i).Score = score;
    alPromSeq(i).Start = startat';
    at = struct;
    at.cer = alignment(1, :);
    at.par = alignment(3, :);
    cerNoGaps = at.cer; cerNoGaps(strfind(cerNoGaps, '-')) = [];
    parNoGaps = at.par; parNoGaps(strfind(parNoGaps, '-')) = [];
    alPromSeq(i).Stop = [startat(1) + length(cerNoGaps)-1, startat(2) + length(parNoGaps)-1];
    alPromSeq(i).Alignment = alignment;
    at_len = length(at.cer);
    for j = 1:2
        gaps = strfind(at.(sp{j}), '-');
        notGaps = 1:at_len; notGaps(gaps) = [];
        currRegion = regionInGenome(j, 1):regionInGenome(j, 2);
        locInAt = alPromSeq(i).Start(j):alPromSeq(i).Stop(j);
%         currRegion = currRegion(locInAt);
        if strand(j) == 0
            currRegion = currRegion(end:-1:1);   
        end
        currRegion = currRegion(locInAt);
        currAlignedLocCer = nan(2, at_len);
        currAlignedLocCer(1, notGaps) = currRegion;
        currAlignedLocCer(2, :) = currGn;
        alignedLoc.(sp{j}){i} = currAlignedLocCer;
    end
    %%
end
%% test
s1 = struct;
currGn = i;
for j = 1:2
    locNotNan = alignedLoc.(sp{j}){currGn}(1, :); 
    locNotNan(isnan(locNotNan)) = [];
    currSeq = wholeGen.(sp{j})(locNotNan);
    if locNotNan(1) > locNotNan(2)
        % minus strand
        currSeq = seqrcomplement(currSeq(end:-1:1));
    end
    s1.(sp{j}) = currSeq;
end

[s, a] = nwalign(s1.cer, s1.par)
%%

save('alignedLoc_1_3.mat', 'alignedLoc');
save('alPromSeq_1_3.mat', 'alPromSeq');
%%

