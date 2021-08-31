
currGP = GP.gene_infoR64.position;
idxNotNan = find(all(~isnan(currGP), 2));
status = GP.gene_infoR64.status;
status(cellfun(@isempty, status)) = {'Dubious ORF'};
DubiousIdx = find(strcmp(status, 'Dubious ORF'));
NotDubious = 1:6701; NotDubious(DubiousIdx) = [];
idxNotNan = intersect(idxNotNan, NotDubious);
currGP = currGP(idxNotNan, :);
strand = nan(length(idxNotNan), 1);
strand(currGP(:, 3) - currGP(:, 2) > 0) = 1;
strand(currGP(:, 3) - currGP(:, 2) < 0) = 0;

sortedGenes = [];
sortedGP = [];
startOfChrs = nan(17, 1);
endsOfChrs = nan(17, 1);
for i = 1:17
    idxGenes = find(currGP(:, 1) == i);
    currChr = currGP(idxGenes, :);
    smallerCoor = min(currChr(:, 2:3), [], 2);
    [~, si] = sort(smallerCoor);
    startOfChrs(i) = idxGenes(si(1));
    endsOfChrs(i) = idxGenes(si(end));
    sortedGenes = [sortedGenes; idxGenes(si)];
    sortedGP = [sortedGP; currChr(si, :)];
end

strand = strand(sortedGenes);
relToPrev = strand - [nan; strand(1:end-1)];
relToPrev(find(ismember(sortedGenes, startOfChrs))) = nan;
head2head = idxNotNan(sortedGenes(find(relToPrev == 1)));
head2head_neighbor = idxNotNan(sortedGenes(find(relToPrev == 1) - 1));
tail2tail = idxNotNan(sortedGenes(find(relToPrev == -1)));
head2tail = idxNotNan(sortedGenes(find(relToPrev == 0)));
disp(['# head2head = ', num2str(length(head2head))]);
disp(['# tail2tail = ', num2str(length(tail2tail))]);
disp(['# head2tail = ', num2str(length(head2tail))]);

% NOTE: head2head means the gene is shares the same intergenic region with
% the upstream gene (for example, GDH3 and FLO9). GDH3 is marked head2head,
% and FLO9 is head2tail with YAL064W.
