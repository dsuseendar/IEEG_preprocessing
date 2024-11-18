function specMean =  applySpecNorm(spec,meanFreqChanIn)

    assert(length(spec)==size(meanFreqChanIn,1),"Channel number mismatch");
    
    for iChan = 1:length(spec)
        spec2Analyze = spec{iChan}; % Spectrogram to analyze   
        meanFreq = meanFreqChanIn(iChan, :); 
        specMean{iChan} = (squeeze(mean(spec2Analyze ./ reshape(meanFreq, 1, 1, []), 1)))'';
       
    end
end