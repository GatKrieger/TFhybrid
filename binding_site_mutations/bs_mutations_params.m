if ~exist('absSumPerPosMutAll', 'var');
    dataDir = [homeDir, 'checSeq_project\analyze\signal_on_aligned_seq\binding_site_mutations\mat_files\'];
    load([dataDir, 'absSumPerPosMutAll-1.mat']);
    load([dataDir, 'LRperPosMutAll-1.mat']);
    load([dataDir, 'seqMatAll-1.mat']);
end

maxLen = max(structfun(@length, PFMtoUse)) + 10;

flankingBases = 5;