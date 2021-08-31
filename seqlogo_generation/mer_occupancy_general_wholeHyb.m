function [out] = mer_occupancy_general_wholeHyb(chec_profile, nmer, varargin)
global genomes tss_struct
% from Felix Jonas
ip=inputParser;
ip.addParameter('window',20);
ip.addParameter('promoterLength',400);
ip.addParameter('downStreamTSS',0);
ip.parse(varargin{:});
window = ip.Results.window;
promoterLength = ip.Results.promoterLength;
downStreamTSS = ip.Results.downStreamTSS;
%%
refGenome = genomes.hyb;
tss = tss_struct.hyb;
intergenicRegion = tss_struct.intergenic.hyb;
Nchr = 32;
disp('generating sc4base');
ntBases={'A','C','G','T'};
sc4base = cell(1, Nchr);
for i=1:Nchr
    [~,sc4base{i}]=ismember(upper(refGenome{i}'),ntBases);
    sc4base{i}=sc4base{i}-1;
end
clear SC_genome
sc4base=cat(1,sc4base{:});
kernel=4.^[0:nmer-1];
disp('finished generating sc4base');
%% create mers table
disp('generating nmer table');
motifVal = cellfun(@str2num, mat2cell(dec2base([0:4^nmer-1]',4, nmer), ones(4^nmer, 1), ones(nmer, 1)));
nmerTable = table(mat2cell(cell2mat(ntBases(motifVal + 1)), ones(4^nmer, 1), nmer), [0:4^nmer-1]' + 1, ...
    'VariableNames',{'seq', 'value'});
% nmerTable.rcValue = sum((3-motifVal).* fliplr(kernel), 2)+1;
nmerTable.rcValue = sum((3-motifVal).* kernel, 2)+1;
nmerTable.rcSeq = nmerTable.seq(nmerTable.rcValue);
sc4nmer = conv(sc4base, kernel', 'same')+1;
%%
nmerRed = nmerTable(nmerTable.value < nmerTable.rcValue, :);
sc4red = changem(sc4nmer, nmerRed.value, nmerRed.rcValue);
sc4red = changem(sc4red,1:size(nmerRed, 1),nmerRed.value);
clear nmerTable sc4bases
disp('finished generating nmer table');
%% Annotating promoters (logical)

% promotersRegion = [];
tmp = []; %temporary promoters regions

for i = 1:Nchr
    tmp{i} = zeros(1,length(refGenome{i}));
end

genes_with_annotated_stanscript = find(~isnan(tss(:,1)));

for i = 1:length(genes_with_annotated_stanscript)
    gene = genes_with_annotated_stanscript(i);
    locations = tss(gene,:);
    strand = sign(locations(3)-locations(2));
    currPromoterLength = intergenicRegion(genes_with_annotated_stanscript(i));
    if isnan(currPromoterLength) | currPromoterLength > promoterLength
        currPromoterLength = promoterLength;
    end
    tmp{locations(1)}((locations(2) - strand*currPromoterLength) : strand : locations(2) + downStreamTSS) = 1;
end

% for i = 1:Nchr
%     promotersRegion = [promotersRegion,tmp{i}];
% end
% promotersRegion = logical(promotersRegion);
promotersRegion = logical([tmp{:}]);
%%
disp('computing scores per profile');
score = zeros(size(nmerRed,1),size(chec_profile,2));

for i=1:size(chec_profile,2)
    disp(num2str(i));
    occ_i = movmean(chec_profile(:, i), 2*window + 1,1);
    score(:, i) = accumarray(sc4red(promotersRegion), occ_i(promotersRegion), [size(nmerRed, 1) 1], @(x)mean(x,'omitnan'));
end
disp('finished');
out = struct;
out.mers = nmerRed;
out.score = score;
out.sc4nmer = sc4nmer;
out.sc4red = sc4red;
out.promotersRegion = promotersRegion;
end