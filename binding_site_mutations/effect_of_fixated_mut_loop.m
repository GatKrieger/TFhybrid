
% factors2 = {'REB1', 'MBP1', 'TBF1', 'RAP1', 'MSN2', 'SUT1', 'MIG1', 'ACE2', 'SWI4', 'SKN7', 'MIG2', 'SWI5', ...
%     'ABF1', 'SOK2', 'PHO4', 'PHD1', 'TOD6', 'COM2', 'FKH1', 'ASH1', 'RGT1', 'DOT6', 'GCR1', 'GCR2', 'SUT2', ...
%     'TEC1', 'YAP2', 'SFP1'};

% factors2 = {'ABF1', 'REB1', 'TOD6', 'RAP1', 'TEC1', 'FKH1', 'SUT1', 'DOT6', 'SWI5', 'MIG1', 'MIG2', ...
%     'MSN2', 'SWI4', 'GCR1', 'SKN7', 'SOK2', 'FKH2', 'RGT1', 'ACE2', 'PHO4', 'GCR2', 'STB3', 'PHD1', 'TBF1', ...
%     'HAP4', 'MBP1', 'HMS2', 'FHL1', 'ASH1', 'PHO2'};

NTF = length(factors2);
tValStruct = struct;
absSumPerPosMutAll = struct;
LRperPosMutAll = struct;
whichSitePerPosMutAll = struct;
rawDataAll = struct;
flankingBases = 5;
mutationRegime = 'onlyOne';
MinPeakHeightByStdfromRandom = 1.5;
for IDX2 = 1:length(factors2)
% for IDX2 = 28:length(factors2)
    currTF = factors2{IDX2};
    pfmMatchScoreTh = minScoreToUse.(currTF);
    disp(currTF);
    tVals = [];
    effect_of_fixated_mut;
    tValStruct.(currTF) = tVals;
    absSumPerPosMutAll.(currTF) = absSumPerPosMut;
    LRperPosMutAll.(currTF) = LRperPosMut;
    whichSitePerPosMutAll.(currTF) = whichSitePerPosMut;
    rawDataAll.(currTF) = rawData;
end

saveDir = 'mat_files_onlyOneMut_inVitro_std1_5/';
save([saveDir, 'absSumPerPosMutAll-1.mat'], 'absSumPerPosMutAll');
save([saveDir, 'LRperPosMutAll-1.mat'], 'LRperPosMutAll');
save([saveDir, 'rawDataAll-1.mat'], 'rawData');

% save('absSumPerPosMutAll-1_oneOrMoreMuts.mat', 'absSumPerPosMutAll');
% save('LRperPosMutAll-1_oneOrMoreMuts.mat', 'LRperPosMutAll');
% save('seqMatAll-1_oneOrMoreMuts.mat', 'seqMatAll');
% save('rawDataAll-1_oneOrMoreMuts.mat', 'rawData');