function pSig = timePermCluster(activeSignal,passiveSignal,clusterargs)
arguments
    activeSignal {mustBeMatrix}
    passiveSignal {mustBeMatrix}
    clusterargs.statstype = 't-stat';
    clusterargs.nPerm = 1000; 
    clusterargs.numTail = 1;
    clusterargs.pThresh = 0.05;
end
pThresh = clusterargs.pThresh;
nPerm = clusterargs.nPerm;
assert(size(activeSignal,2)==size(passiveSignal,2),'Signals have different time length');

statMetric = statfun(activeSignal,passiveSignal,type = clusterargs.statstype);

combSignals = cat(1,activeSignal,passiveSignal);
nComb = size(combSignals,1);
fHalf = 1:size(activeSignal,1);
sHalf = size(activeSignal,1)+1:nComb;
fprintf(1, '\n>> Calculating permutation statistic...\n');
for iPerm = 1:nPerm
    permIds = randperm(nComb);
    permStat(iPerm,:) = statfun(combSignals(permIds(fHalf),:),combSignals(permIds(sHalf),:),type = clusterargs.statstype);
    fprintf(1,'.');
end
fprintf(1,'] done\n')

pStatOne = sum(permStat > statMetric,1)./nPerm;
pStatTwo = sum(abs(permStat) > abs(statMetric),1)./nPerm;

pStatOne = adjustPVals(pStatOne,nPerm);
pStatTwo = adjustPVals(pStatTwo,nPerm);

fprintf(1, '\n>> Calculating permutation p-values...\n');
for iPerm = 1:nPerm
    permRecur = setdiff(1:nPerm, iPerm);
    pStatPermOne(iPerm,:) = sum(permStat(iPerm,:)>permStat(permRecur,:),1)./(nPerm-1);
    pStatPermOne(iPerm,:) = adjustPVals(pStatPermOne(iPerm,:),nPerm-1);
    pStatPermTwo(iPerm,:) = sum(abs(permStat(iPerm,:))>abs(permStat(permRecur,:)),1)./(nPerm-1);
    pStatPermTwo(iPerm,:) = adjustPVals(pStatPermTwo(iPerm,:),nPerm-1);
    fprintf(1,'.');
end
fprintf(1,'] done\n')
if(clusterargs.numTail == 1)
    fprintf(1,'\n>> Performing one-sided cluster correction')
    [pValsRaw, actClust] = timePermClusterAfterPermPValues(pStatOne, pStatPermOne, pThresh);
else
    fprintf(1,'\n>> Performing two-sided cluster correction')
    [pValsRaw, actClust] = timePermClusterAfterPermPValues(pStatTwo, pStatPermTwo, pThresh);

end

h_sig_05 = zeros(1,length(pValsRaw));
h_sig_01 = zeros(1,length(pValsRaw));
for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        h_sig_05(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
        if actClust.Size{iClust}>actClust.perm99
            h_sig_01(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end
end

pSig.pVals = pValsRaw;
pSig.clust = actClust;
pSig.h_sig_05 = h_sig_05;
pSig.h_sig_01 = h_sig_01;

end

function statMetric = statfun(activeSignal, passiveSignal, statargs)
arguments
    activeSignal % trials x time
    passiveSignal % trials x time
    statargs.type {mustBeMember(statargs.type, {'mean-sub', 't-stat', 'f-stat'})} = 'mean-sub'
end

switch statargs.type
    case 'mean-sub'
        statMetric = mean(activeSignal, 1) - mean(passiveSignal, 1);
    
    case 't-stat'
        [~, ~, ~, stats] = ttest2(activeSignal, passiveSignal);
        statMetric = stats.tstat;
    
    case 'f-stat'
        % Compute F-statistic
        n1 = size(activeSignal, 1);
        n2 = size(passiveSignal, 1);
        var1 = var(activeSignal);
        var2 = var(passiveSignal);
        statMetric = var1 ./ var2;
        
        % Adjust F-statistic if var1 < var2
        flip_idx = var1 < var2;
        statMetric(flip_idx) = 1 ./ statMetric(flip_idx);
end

end

function [adjustedPVals] = adjustPVals(pVals, nPerm)
    adjustedPVals = max(min(pVals, 1 - 1/nPerm), 1/nPerm);
    
end