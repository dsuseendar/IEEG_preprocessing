function ieegPad = extendTimeEpoch(ieegData,sigLen)
%EXTENDTIMEEPOCH 
warning off;
padTimeArray = 20:20:size(ieegData,3);


for iChan = 1:size(ieegData,1)
    ieegChan = squeeze(ieegData(iChan, :, :));
    ieegChanPad = [];
    for iTrial = 1:size(ieegChan, 1)
        padChoose = randsample(padTimeArray,1);
        time2pad = sigLen/padChoose;
        selectTrials = setdiff(1:size(ieegChan, 1),iTrial);
        randTrials = datasample(selectTrials, ceil(time2pad) - 1, 'Replace', false);

        trials2join = ieegChan(randTrials, 1:padChoose)';
        %sigGenLen = length(trials2join(:)) +padChoose;
        if(ceil(time2pad)==time2pad)
            ieegChanPad(iTrial, :) = [ieegChan(iTrial, 1:padChoose) trials2join(:)'];
        else
            %time2remove = (ceil(time2pad)-(time2pad))*fs;
            joinTrials = trials2join(:)';
            %randTrial = ieegChan(randsample(selectTrials,1),1:extraTimeNeed*fs-1);
            ieegChanPad(iTrial, :) = [ieegChan(iTrial, 1:padChoose) joinTrials(1:sigLen-padChoose)];
        end
    end
%     for iTrial = 1:size(ieegChan, 1)
%         %ieegChanPad(iTrial, :) = ieegChan(iTrial, :);
%         randStart = randi(20);
%         ieegflip = ieegChan(iTrial, randStart:randStart+padNum-1);
%         ieegPadTrial= ieegChan(iTrial, randStart:randStart+padNum-1);
%         for iPad = 1:time2pad-1 
%             ieegflip = fliplr(ieegflip);
%             ieegPadTrial = [ieegPadTrial; ieegflip];
%         end
%         ieegChanPad(iTrial,:) = ieegPadTrial;
%     end
    ieegPad(iChan,:,:) = ieegChanPad;
end



