disp('score motif by Kmers');
KmerScoreS = struct;
alMotifScoresS = struct;
motifTechnique = 2;
scoreBy = 'av';
Xmer = 7;
for i = 1:length(factors)
    currTF = factors{i};
    disp(currTF);
    score_motif_by_kmers;
    KmerScoreS.(currTF) = KmerScore;
    alMotifScoresS.(currTF) = alMotifScores;
end

save('alMotifScoresS.mat', 'alMotifScoresS');
save('KmerScoreS.mat', 'KmerScoreS');