function [spec, pPerc] = getSpectrograms(ieeg, goodtrials, tw, etw, efw, prtw, pertw, intF, fs, ispermTest)
% getSpectrograms - Extracts spectrograms and performs statistical tests on ECoG data.
%
% Inputs:
%    ieeg - ECoG data (channels x trials x time)
%    goodtrials - Good trials for each channel (cell array of size [1 x channels])
%    tw - Time window of interest [start_time, end_time] in seconds
%    etw - Spectrogram time window [start_time, end_time] in seconds
%    efw - Spectrogram frequency window [start_frequency, stop_frequency] in Hz
%    prtw - Pre-onset time window [start_time, end_time] to get significant channels
%    pertw - Post-onset time window [start_time, end_time] to get significant channels
%    intF - Frequency range of interest for statistical tests [start_frequency, stop_frequency] in Hz
%    fs - Sampling frequency in Hz
%    ispermTest - Flag (0/1) indicating whether to perform a permutation test to determine channel significance
%
% Outputs:
%    spec - Spectrograms of each trial for each channel (cell array of size [1 x channels])
%    pPerc - P-values from the permutation test to check channel significance (1 x channels)
%
% Example:
%    ieeg = rand(10, 100, 1000); % Example ECoG data
%    goodtrials = cell(1, 10); % Example good trials (cell array)
%    tw = [0, 10]; % Time window of interest
%    etw = [2, 8]; % Spectrogram time window
%    efw = [30, 80]; % Spectrogram frequency window
%    prtw = [0, 2]; % Pre-onset time window
%    pertw = [2, 4]; % Post-onset time window
%    intF = [30, 50]; % Frequency range of interest for statistical tests
%    fs = 1000; % Sampling frequency
%    ispermTest = 1; % Perform permutation test
%    [spec, pPerc] = getSpectrograms(ieeg, goodtrials, tw, etw, efw, prtw, pertw, intF, fs, ispermTest); % Extract spectrograms and perform statistical tests
%

 % Parameters for spectrogram analysis
 AnaParams.dn = 0.05;
 AnaParams.Tapers = [0.5, 10];
 AnaParams.fk = [efw(1), efw(2)];
 AnaParams.Fs = fs;

 % Channels and permutation parameters
 channelOfInterest = 1:size(ieeg, 1);
 numPerm = 10000;

 % Arrayfun to process each channel and store results
 [spec, pPerc] = arrayfun(@(iChan) process_channel(iChan, ieeg, goodtrials, AnaParams, ...
                                                   tw, etw, prtw, pertw, intF, ispermTest, numPerm), ...
                          channelOfInterest, 'UniformOutput', false);

 % Convert pPerc from cell array to numeric array for output
 pPerc = cell2mat(pPerc);
end

% Helper function to process each channel
function [spec_out, pPerc_out] = process_channel(iChan, ieeg, goodtrials, AnaParams, ...
                                              tw, etw, prtw, pertw, intF, ispermTest, numPerm)
 % Determine trials for current channel
 if isempty(goodtrials)
     trials_g = 1:size(ieeg, 2);
 elseif iscell(goodtrials)
     trials_g = goodtrials{iChan};
 else
     trials_g = goodtrials;
 end

 % Extract spectrogram and define frequency/time windows
 [spec_out, F] = extract_spectrograms_channel(squeeze(ieeg(iChan, trials_g, :)), AnaParams);
 gammaFreq = F >= intF(1) & F <= intF(2);
 tspec = linspace(tw(1), tw(2), size(spec_out, 2));
 prtspec = tspec >= prtw(1) & tspec <= prtw(2);
 perctspec = tspec >= pertw(1) & tspec <= pertw(2);

 % Perform permutation test if required
 if ispermTest == 1
     meanBase = arrayfun(@(t) mean2(squeeze(spec_out(t, prtspec, gammaFreq))), 1:length(trials_g));
     meanOnsetPercept = arrayfun(@(t) mean2(squeeze(spec_out(t, perctspec, gammaFreq))), 1:length(trials_g));
     pPerc_out = permtest(meanOnsetPercept, meanBase, numPerm);
 else
     pPerc_out = 0;
 end

 % Extract time window of interest
 etspec = tspec >= etw(1) & tspec <= etw(2);
 spec_out = spec_out(:, etspec, :);
end
