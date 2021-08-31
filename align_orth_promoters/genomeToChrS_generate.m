%% genome to chromosomes
genomeToChrS = struct;
for j = 1:2
    currGenome = genomes.(spLongName{j});
    chrLen = cellfun(@length, currGenome);
    genomeLen = sum(chrLen);
    genomeToChr = cell(length(chrLen), 1);
    for i = 1:length(chrLen)
        genomeToChr{i} = [1:chrLen(i); i*ones(1, chrLen(i))];
    end
    genomeToChr = [genomeToChr{:}];
    genomeToChrS.(sp{j}) = genomeToChr;
end
save('genomeToChrS.mat', 'genomeToChrS');