averageRepeats = 1;
totalbasesAroundPeak = 61;
basesAroundPeak = (totalbasesAroundPeak - 1)/2;
basesAroundPeakForSum = 10;
prctileToFlt = 90;
confRange = 100;
ifScale = 1;
scaleBy = 'sumOnPromoter';
% scaleBy = 'fullData';
scoreBy = 'av';
Xmer = 7;
motifTechnique = 2;
MinPeakDistance = 20;
% MinPeakDistance = 50;
minNreads = 5;
smoothingSize = 20;
MinPeakWidth = 10;
Kpeaks = 1000;
restrictKpeaks = 0;
LRth = 1;
distFromTSSTh = 800;
% MinPeakHeightFixed = 500;
L = length(wholeAlignment);
currSample=currTF;
if ~exist('tss_wholeGenIdx', 'var')
    load([homeDir, 'checSeq_project/analyze/signal_on_aligned_seq/data_structs_swalign/tss_wholeGenIdx.mat']);
end