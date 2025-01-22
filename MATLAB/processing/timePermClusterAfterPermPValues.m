function [pValsRaw, actClust] = timePermClusterAfterPermPValues(testP, permP, pThresh)
% timePermClusterAfterPermPValues - Compute time-based permutation clustering analysis.
% 
% Inputs:
%   testP: P-values for the test signal (1 x time).
%   permP: P-values for permuted signals (permutations x time).
%   pThresh: Threshold for significance (e.g., 0.05).
% 
% Outputs:
%   pValsRaw: Raw p-values for the actual data.
%   actClust: Cluster information.

% Validate input dimensions
if size(testP, 2) ~= size(permP, 2)
    error('Signals have different time lengths');
end

% Get the number of permutations
nPerm = size(permP, 1);

% Compute clusters for actual data
tmpA = testP < pThresh;
cc = bwconncomp(tmpA);

% Initialize cluster information structure
%actClust = struct('Start', {}, 'Size', {}, 'maxPermClust', [], 'perm95', [], 'perm99', []);

if ~isempty(cc.PixelIdxList)
    for ii = 1:numel(cc.PixelIdxList)
        actClust.Start{ii} = cc.PixelIdxList{ii}(1);
        actClust.Size{ii} = numel(cc.PixelIdxList{ii});
    end
else
    actClust.Start = {NaN};
    actClust.Size = {NaN};
end

% Compute maximum cluster sizes for permuted data
maxPermClust = zeros(nPerm, 1);
fprintf(1, '\n>> Identifying connected time-points in permuted values\n');
fprintf(1,'[');
parfor iPerm = 1:nPerm
    tmpA = permP(iPerm, :) < pThresh;
    cc = bwconncomp(tmpA);
    if ~isempty(cc.PixelIdxList)
        maxPermClust(iPerm) = max(cellfun(@numel, cc.PixelIdxList));
    end
    fprintf(1,'.');
end
fprintf(1,'] done\n');

% Calculate percentiles
actClust.maxPermClust = sort(maxPermClust);
actClust.perm95 = prctile(actClust.maxPermClust, 95);
actClust.perm99 = prctile(actClust.maxPermClust, 99);

% Return raw p-values
pValsRaw = testP;
end
