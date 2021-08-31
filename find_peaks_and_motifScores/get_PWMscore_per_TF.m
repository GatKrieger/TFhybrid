% get PWMscore for a TF
cerGenLen = 12071326;
hybGenLen = 24092654;
genInHyb = {[1, cerGenLen], [cerGenLen+1 hybGenLen]};
pfm = PFMtoUse.(currTF);
motif_length = size(pfm, 2);
if ~strcmp(currTF, 'ABF1')
    nmerStr = ['mer', num2str(motif_length)];
    pfmScore = pbsPwm.(currTF).pfmScore;
    currAlGenXmer = struct;
    alPWMscore = struct;
    genPWMscore = struct;
    Nopts = 4^motif_length;
    for innerIdx = 1:2
        currSpi_nnIdx = sp{innerIdx};
        currAlGenXmer.(currSpi_nnIdx) = alGenomeXmers.(nmerStr).(currSpi_nnIdx);
        curr_alPWMscore = nan(1, L);
        idxNotNan = find(~isnan(currAlGenXmer.(currSpi_nnIdx)));
        curr_alPWMscore(idxNotNan) = pfmScore(currAlGenXmer.(currSpi_nnIdx)(idxNotNan));
        alPWMscore.(currSpi_nnIdx) = curr_alPWMscore;
        
        currGenXmer = KmerTables.(nmerStr).genomeXmers(genInHyb{innerIdx}(1) : genInHyb{innerIdx}(2));
        currGenXmer(currGenXmer > Nopts) = Nopts;
        currGenXmer(currGenXmer < 1) = 1;
        currGenPWMScore = pfmScore(currGenXmer);
        genPWMscore.(currSpi_nnIdx) = currGenPWMScore;
    end
    
else
    genPWMscore = struct;
    genPWMscore.cer = pbsPwm.ABF1.pbsScore_onGen(1:cerGenLen);
    genPWMscore.par = pbsPwm.ABF1.pbsScore_onGen(cerGenLen+1:end);
    alPWMscore = struct;
    
    for innerIdx = 1:2
        currSpi_nnIdx = sp{innerIdx};
        idxNotNan = find(~isnan(alignedLocV.(currSpi_nnIdx)));
        curr_alPWMscore = nan(1, L);
        curr_alPWMscore(idxNotNan) = genPWMscore.(currSpi_nnIdx)(idxNotNan);
        alPWMscore.(currSpi_nnIdx) = curr_alPWMscore;
    end
end
