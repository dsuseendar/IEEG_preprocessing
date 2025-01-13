function meanFreqChanOut = extractSpecNorm(spec, tw, etw)
% extractSpecNorm - Extracts the normalized mean frequency from spectrograms.
%
% Inputs:
%    spec - Spectrograms (cell array of size [1 x nChannels])
%    tw - Time window of interest [start_time, end_time]
%    etw - Extraction time window [start_time, end_time]
%
% Output:
%    meanFreqChanOut - Normalized mean frequency (nChannels x nFrequencies)
%
% Example:
%    spec = cell(10, 1); % Example spectrograms (cell array)
%    tw = [0, 10]; % Time window of interest
%    etw = [2, 8]; % Extraction time window
%    meanFreqChanOut = extractSpecNorm(spec, tw, etw); % Extract normalized mean frequency
%


tspec = linspace(tw(1), tw(2), size(spec{1}, 2)); % Time vector for the time window of interest


% Use arrayfun to calculate mean frequencies for each channel and frequency
meanFreqChanOut = arrayfun(@(iChan) ...
    arrayfun(@(iFreq) ...
        mean2(squeeze(spec{iChan}(:, tspec >= etw(1) & tspec <= etw(2), iFreq))), ...
    1:size(spec{iChan}, 3)), ...
1:length(spec), 'UniformOutput', false);

% Convert the cell array to a matrix
meanFreqChanOut = cell2mat(meanFreqChanOut');

end
