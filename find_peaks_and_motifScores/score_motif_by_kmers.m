
currMer = ['mer', num2str(Xmer)];
KmerScore = motifStruct.hyb.mer7.meanSignalPerMotif.(currTF);

KmerScore2 = nan(length(KmerScore)*2, size(KmerScore, 2));
KmerScore2(motifStruct.hyb.(currMer).table.value, :) = KmerScore;
KmerScore2(motifStruct.hyb.(currMer).table.rcValue, :) = KmerScore;
KmerScore = KmerScore2;
clear KmerScore2;

alMotifScores = struct;
for i = 1:length(sp)
    currSp = sp{i};
    msVec = nan(1, L);
    for j = 1:length(KmerScore)
        msVec(alGenomeXmers.(currSp) == j) = KmerScore(j);
    end
    alMotifScores.(currSp) = msVec;
end

