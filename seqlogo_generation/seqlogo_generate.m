function [finalFreq,totAGCT,Positions2Plot]= seqlogo_generate(motData, MotSet, strain, MotifSeqs)
% from Sagie Brodsky
meanSignal_sample =  motData.(strain);
strainName = strrep(strain, '_', ' ');
motNum = length(MotSet);

MotStet_idx=[];
for m=1:length(MotSet)
    m_idx = find(strcmp({MotSet{m}},MotifSeqs));
    MotStet_idx=[MotStet_idx,m_idx];
end

[meanSignal_sorted, meanSignal_sorted_idx] =sort(meanSignal_sample(MotStet_idx), 'descend');
meanSignal_sorted=meanSignal_sorted(1:motNum);
meanSignal_sorted_idx=meanSignal_sorted_idx(1:motNum);
MotSet = MotSet(meanSignal_sorted_idx);
meanSignal_sorted  = round(meanSignal_sorted);
motRank  = [];
for i = 1:length(MotSet)
    motRank(i) = round((meanSignal_sorted(i) / sum(meanSignal_sorted)),2) * 100;
end

motRank = int64(motRank);

for i = 1:length(MotSet)
    forAlign = nwalign(MotSet{1},MotSet{i},'glocal',true);
    revAlign = nwalign(MotSet{1},seqrcomplement(MotSet{i}),'glocal',true);
    if revAlign > forAlign
        MotSet{i} = seqrcomplement(MotSet{i});
    end
end

MotSetOut= MotSet;

finalMatforAlign = cell(sum(motRank)*2,1);
j=0;
for i = 1:2:sum(motRank )*2-1
    j=j+1;
    finalMatforAlign{i} = ['>',num2str(j)];
end
e = 0;
for f = 1:length(motRank )
    currSize = motRank (f);
    for i = 1:currSize
        finalMatforAlign{int64((e+i)*2)} = MotSet{f};
    end
    e = e + motRank (f);
end

formatAlign = fastaread(char(finalMatforAlign));
finalAlign = multialign(formatAlign,'terminalGapAdjust',true);
sumA = [];
sumG = [];
sumC = [];
sumT = [];
totAGCT = [];
for h= 1:length(finalAlign(1).Sequence)
    currPos = char();
    for i = 1:length(finalAlign)
        currPos = [currPos,char(finalAlign(i).Sequence(h))];
        sumA(h) = sum(currPos=='A');
        sumG(h) = sum(currPos=='C');
        sumC(h) = sum(currPos=='G');
        sumT(h) = sum(currPos=='T');
        totAGCT(h) = length(regexp(currPos,'[AGCT]'));
    end
end

sumAfreq = sumA ./ totAGCT;
sumGfreq = sumG ./ totAGCT;
sumCfreq = sumC ./ totAGCT;
sumTfreq = sumT ./ totAGCT;
finalFreq = [sumAfreq;sumGfreq;sumCfreq;sumTfreq];

[~,iMax] = max(smooth(totAGCT));
Positions2Plot = [iMax-4:iMax+4];

% mySeqLogo(finalFreq(:,Positions2Plot(1):Positions2Plot(end)))
%%
end
