Pos = fields(LRperPosMut);
Npos = length(Pos);
motif_length = length(origConsSeq);
motifLoc = [flankingBases+1, flankingBases+motif_length];
motifLims = [motifLoc(1)-.5 motifLoc(2)+.5 ];
deltaBinding = nan(4, Npos, 2);
Freq = nan(4, Npos);
% minSites = 3;
for i = 1:Npos
    if strcmp(toPlot, 'abs')
%         currData = absSumPerPosMut.(Pos{i}).inGen;
        currData = absSumPerPosMut.(Pos{i}).inGenConsLR;
        Freq(:, i) = cellfun(@length, absSumPerPosMut.(Pos{i}).inGen(2:end));
    elseif strcmp(toPlot, 'LR')
        currData = LRperPosMut.(Pos{i}).bwSp;
        Freq(:, i) = cellfun(@length, currData(2:end));
    end
    currData = currData(2:end); % first is INDEL, not relevant here
    currData(Freq(:, i) < minSites) = {nan};
    deltaBinding(:, i, 1) = cellfun(@nanmedian, currData);
    deltaBinding(:, i, 2) = cellfun(@nanstd, currData);
end
colors = [rgb('Green'); rgb('Blue'); 0.9290    0.6940    0.1250; rgb('Red')];

% figure; 
hold on;
for i = 1:4
    scatter(1:Npos, deltaBinding(i, :, 1), [], colors(i, :), 'filled');
%     errorbar(1:Npos, d(i, :, 1), d(i, :, 2), '.', 'color', colors(i, :));
end
xlim([1, Npos]);
ylim([-6 1]);
axis square; set(gcf, 'color', 'w');
plot([motifLims; motifLims], ylim, '--k');
ylabel({'Change in binding:',  'log2 (alternative / consensus)'});
xlabel('Position');
xlabels = horzcat(cellfun(@num2str, num2cell(-flankingBases:-1), 'uniformoutput', false), ...
    cellfun(@num2str, num2cell(1:flankingBases+motif_length), 'uniformoutput', false));
set(gca, 'xtick', 1:Npos, 'xticklabels', xlabels);
legend({'A', 'C', 'G', 'T'});
for i = 1:motif_length
    text(motifLoc(1) + i-1, 1, origConsSeq(i));
end
grid on;