function [signalPerChr] = bedGraph2out(currBedGraphFile, chrLengths)

bgData = importdata(currBedGraphFile);
chrs = categorical(bgData.textdata);
chrNames = unique(bgData.textdata, 'stable');
Nchr = length(chrLengths);
signalPerChr = cell(Nchr, 1);
for i = 1:Nchr
    idx = find(chrs == chrNames{i});
    currChr = zeros(1, chrLengths(i));
    for j = 1:length(idx)
        currPos = bgData.data(idx(j), 1)+1 : bgData.data(idx(j), 2);
        currChr(currPos) = repmat(bgData.data(idx(j), 3), length(currPos), 1);
    end
    currChr = currChr(1:chrLengths(i));
    signalPerChr{i} = currChr;
end
end
    