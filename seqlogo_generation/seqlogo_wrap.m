function [finalFreq,totAGCT,Positions2Plot] = seqlogo_wrap(checStruct, strain, motNum, currMer, varargin)
% from Sagie Brodsky
p = inputParser();
p.addParameter('subplotDirections', []);
p.parse(varargin{:});
spDir = p.Results.subplotDirections;

motData = checStruct.(currMer).meanSignalPerMotif;
AllMotifs_7mer = checStruct.(currMer).motifs_seq;
[~, topIdx] = maxk(motData.(strain), motNum);
MotSet = AllMotifs_7mer(topIdx);
[finalFreq, totAGCT,Positions2Plot]= seqlogo_Sagie(motData, MotSet, strain, AllMotifs_7mer);

pfm = finalFreq(:,Positions2Plot(1):Positions2Plot(end));
mySeqLogo(pfm, 'subplotDirections', spDir);
end
