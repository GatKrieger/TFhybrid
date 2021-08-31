% get PWMscore for a TF: similar to "get_PWMscore_per_TF" but here using
% the PWM defined by our data (hyb)
cerGenLen = 12071326;
hybGenLen = 24092654;
genInHyb = {[1, cerGenLen], [cerGenLen+1 hybGenLen]};
pfm = PFMhyb.(currTF);
motif_length = size(pfm, 2);
nmerStr = ['mer', num2str(motif_length)];
pfmScore = pbsPwmHyb.(currTF).pfmScore;
currAlGenXmer = struct;
alPWMscore_hyb = struct;
genPWMscore_hyb = struct;
Nopts = 4^motif_length;
for innerIdx = 1:2
    currSpi_nnIdx = sp{innerIdx};
    currAlGenXmer.(currSpi_nnIdx) = alGenomeXmers.(nmerStr).(currSpi_nnIdx);
    curr_alPWMscore = nan(1, L);
    idxNotNan = find(~isnan(currAlGenXmer.(currSpi_nnIdx)));
    curr_alPWMscore(idxNotNan) = pfmScore(currAlGenXmer.(currSpi_nnIdx)(idxNotNan));
    alPWMscore_hyb.(currSpi_nnIdx) = curr_alPWMscore;
    
    currGenXmer = KmerTables.(nmerStr).genomeXmers(genInHyb{innerIdx}(1) : genInHyb{innerIdx}(2));
    currGenXmer(currGenXmer > Nopts) = Nopts;
    currGenXmer(currGenXmer < 1) = 1;
    currGenPWMScore = pfmScore(currGenXmer);
    genPWMscore_hyb.(currSpi_nnIdx) = currGenPWMScore;
end

