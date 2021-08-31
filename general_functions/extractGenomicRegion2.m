function [data, loc] = extractGenomicRegion2(gn, wantedRegion, upStreamBases, downStreamBases, currGP, currData, ifFigure) 
% downStreamBases = downStreamBases-1;
%%
pos = currGP(gn, :);
if all(isnan(pos))
    error('no position');
elseif pos(1)==17
    error('mitochondrial genome');
end
currChrLen = length(currData{pos(1)});
ORFlength = abs(pos(3) - pos(2));
if pos(2) < pos(3)
    strand = '+';
    % end of Chromosome
    if pos(2) < upStreamBases
        upStreamBases = pos(2)-1;
    end
    if pos(3)+downStreamBases > currChrLen
        downStreamBases = currChrLen - pos(3)-1;
    end
    
    if strcmp(wantedRegion, 'start')
        region = [pos(2) - upStreamBases, pos(2) + downStreamBases];
        x = -upStreamBases:downStreamBases;
    elseif strcmp(wantedRegion, 'orf')
        region = [pos(2) - upStreamBases, pos(3) + downStreamBases];
        x = -upStreamBases:ORFlength + downStreamBases;
    elseif strcmp(wantedRegion, 'end')
        region = [pos(3) - upStreamBases, pos(3) + downStreamBases];
        x = -upStreamBases:downStreamBases;
    end
    data = currData{pos(1)}(region(1):region(2));
else
    strand = '-';
    % end of Chromosome
    if pos(2) < downStreamBases
        downStreamBases = pos(2)-1;
    end
    if pos(2)+upStreamBases > currChrLen
        upStreamBases = currChrLen - pos(2)-1;
    end
    
    if strcmp(wantedRegion, 'start')
        region = [pos(2)- downStreamBases, pos(2) + upStreamBases];
        x = -upStreamBases:downStreamBases;
    elseif strcmp(wantedRegion, 'orf')
        region = [pos(3)- downStreamBases, pos(2) + upStreamBases ];
        x = -upStreamBases:ORFlength + downStreamBases;
    elseif strcmp(wantedRegion, 'end')
        region = [pos(3)- downStreamBases, pos(3) + upStreamBases];
        x = -upStreamBases:downStreamBases;
    end
    data = currData{pos(1)}(region(1):region(2));
    if ischar(data)
        data = seqrcomplement(data);
    elseif isnumeric(data)
        data = data(end:-1:1);
    end
end
if ifFigure
    figure; plot(x, data)
    xlim([x(1) x(end)]);
end
if size(data, 1) ~= 1
    data = data';
end
loc = [pos(1) region];
%%
end