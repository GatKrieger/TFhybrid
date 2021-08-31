%% signal around motif
signalAtBS = struct;
signalAtBS_linear = struct;
NucAtBS = struct;
basesAroundPeakForPlot = 100;
totalbasesAroundPeakForPlot = 2*basesAroundPeakForPlot;
rawFolder = 'signal_around_motif_raw/';
for IDX = 1:NTF
    currTF = reorderedfactors2{IDX};
    currCerPBS = [pbsPwm.(currTF).loc_cer_plus_alG; pbsPwm.(currTF).loc_cer_minus_alG]; 
    currParPBS = [pbsPwm.(currTF).loc_par_plus_alG; pbsPwm.(currTF).loc_par_minus_alG]; 
    potentialBS = [currCerPBS; currParPBS];
    whichSp = [ones(length(currCerPBS), 1); 2*ones(length(currParPBS), 1)];
    pbsStrand = [ones(length(pbsPwm.(currTF).loc_cer_plus_alG), 1); ...
        zeros(length(pbsPwm.(currTF).loc_cer_minus_alG), 1); ...
        ones(length(pbsPwm.(currTF).loc_par_plus_alG), 1); ...
        zeros(length(pbsPwm.(currTF).loc_par_minus_alG), 1)];
    NBS = length(potentialBS);
    idxTF = find(strcmp(factors, currTF));
    currSigSm = [signalAlignedSm.cer(:, idxTF), signalAlignedSm.par(:, idxTF)];
    currSig = [signalAligned.cer(:, idxTF), signalAligned.par(:, idxTF)];

    get_signal_at_potential_bs;
    
    currsigmasIdx_1 = sigmasIdx_1(si1_5(IDX), :);
    currsigmasIdx_2 = sigmasIdx_2(si1_5(IDX), :);
    currsigmasIdx_1(isnan(currsigmasIdx_1)) = [];
    currsigmasIdx_2(isnan(currsigmasIdx_2)) = [];
    
    signalAroundMotifB = [signalAroundMotifA(whichSp==1, :, 1); ...
        signalAroundMotifA(whichSp==2, :, 2)];
    signalAroundMotifBLinear = [signalAroundMotifALinear(whichSp==1, :, 1); ...
        signalAroundMotifALinear(whichSp==2, :, 2)];
    NucAroundMotifB = [NucAroundMotifA(whichSp==1, :, 1); ...
        NucAroundMotifA(whichSp==2, :, 2)];
    signalAroundMotifB(pbsStrand == 0, :) = fliplr(signalAroundMotifB(pbsStrand == 0, :));
    signalAroundMotifBLinear(pbsStrand == 0, :) = fliplr(signalAroundMotifBLinear(pbsStrand == 0, :));
    NucAroundMotifB(pbsStrand == 0, :) = fliplr(NucAroundMotifB(pbsStrand == 0, :));
    
    signalAroundMotifB = signalAroundMotifB(sortedSitesIdx.(currTF), :);
    signalAroundMotifBLinear = signalAroundMotifBLinear(sortedSitesIdx.(currTF), :);
    NucAroundMotifB = NucAroundMotifB(sortedSitesIdx.(currTF), :);
    
    % save raw data
    signalAM = signalAroundMotifBLinear;
    save([rawFolder, 'signalAM', currTF, '.mat'], 'signalAM');
    potentialBS = randomSites.(currTF);
    get_signal_at_potential_bs;
    signalAMRandom = [signalAroundMotifALinear(whichSp==1, :, 1); ...
        signalAroundMotifALinear(whichSp==2, :, 2)];
    save([rawFolder, 'signalAMRandom', currTF, '.mat'], 'signalAMRandom');
    save([rawFolder, 'NucAroundMotifB' currTF, '.mat'], 'NucAroundMotifB');
   
    m = nan(Nsigmas, totalbasesAroundPeakForPlot);
    mLinear = m;
    mNuc = nan(Nsigmas, totalbasesAroundPeakForNuc);
    for J = 1:length(currsigmasIdx_1)
        currIdx = currsigmasIdx_1(J) : currsigmasIdx_2(J);
        if length(currIdx) < 5; continue; end
        m(J, :) = nanmean(signalAroundMotifB(currIdx, :), 1);
        mLinear(J, :) = nanmean(signalAroundMotifBLinear(currIdx, :), 1);
        mNuc(J, :) = nanmean(NucAroundMotifB(currIdx, :), 1);
    end
    signalAtBS.(currTF) = m;
    signalAtBS_linear.(currTF) = mLinear;
    NucAtBS.(currTF) = mNuc;
end

%% plot peak around motif
ax4 = nexttile([1 2]);
xlabels = -basesAroundPeakForPlot:basesAroundPeakForPlot:basesAroundPeakForPlot;
xlabels = cellfun(@num2str, num2cell(xlabels), 'uniformoutput', false);
xlocs = [1, basesAroundPeakForPlot, totalbasesAroundPeakForPlot];
m = cell2mat(struct2cell(signalAtBS));
imagesc(m); 
set(gca, 'xtick', xlocs, 'xticklabels', xlabels);
hold on; plot(xlim, [Ncats + .5:Ncats:NTF*Ncats; Ncats + .5:Ncats:NTF*Ncats]', 'w');
set(gca, 'ytick', Ncats-1:Ncats:NTF*Ncats, 'yticklabels', reorderedfactors2);
cb = colorbar('southoutside'); cb.Label.String = 'log2 ChEC-seq signal'; cb.FontSize = fs;
currCm = cbrewer2('BuPu'); 
colormap(ax4, currCm)
caxis([0 4]);
rectPos = nan(NTF, 4);
for i = 1:NTF
    currTF = reorderedfactors2{i};
    motif_length = size(PFMtoUse.(currTF), 2);
    halfMotifLen1 = floor(motif_length ./ 2);
    halfMotifLen2 = motif_length - halfMotifLen1;
    toAdd=1;
    if rem(motif_length, 2) == 0
        toAdd = 0;
    end
    basesBeforePeakPlus = basesAroundPeakForPlot - halfMotifLen1 + toAdd;
    rectPos(i, :) = [basesBeforePeakPlus-.5, .5+i*3-3, motif_length, 3];
    rectangle('Position', rectPos(i, :));
end

%% plot nucleosomes around motif
ax5 = nexttile([1 2]);
xlabels = -basesAroundPeakForNuc:basesAroundPeakForNuc:basesAroundPeakForNuc;
xlabels = cellfun(@num2str, num2cell(xlabels), 'uniformoutput', false);
xlocs = [1, basesAroundPeakForNuc, totalbasesAroundPeakForNuc];
m = cell2mat(struct2cell(NucAtBS));
imagesc(m);
set(gca, 'xtick', xlocs, 'xticklabels', xlabels);
hold on; plot(xlim, [Ncats + .5:Ncats:NTF*Ncats; Ncats + .5:Ncats:NTF*Ncats]', 'w');
set(gca, 'ytick', Ncats-1:Ncats:NTF*Ncats, 'yticklabels', reorderedfactors2);
cb = colorbar('southoutside'); cb.Label.String = 'Nucleosome occupancy'; cb.FontSize = fs;
currCm = cbrewer2('GnBu'); 
colormap(ax5, currCm)
cb = colorbar('southoutside'); cb.Label.String = 'Nucleosome occupancy'; cb.FontSize = fs;
for i = 1:NTF
    currTF = reorderedfactors2{i};
    motif_length = size(PFMtoUse.(currTF), 2);
    halfMotifLen1 = floor(motif_length ./ 2);
    halfMotifLen2 = motif_length - halfMotifLen1;
    toAdd=1;
    if rem(motif_length, 2) == 0
        toAdd = 0;
    end
    basesBeforePeakPlus = basesAroundPeakForNuc - halfMotifLen1 + toAdd;
    rectPos(i, :) = [basesBeforePeakPlus-.5, .5+i*3-3, motif_length, 3];
    rectangle('Position', rectPos(i, :));
end