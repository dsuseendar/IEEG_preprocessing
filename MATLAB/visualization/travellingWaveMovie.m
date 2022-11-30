function [waveChannel,timeSelect] = travellingWaveMovie(sig2Movie,chanMap,timeAll,etw,clim,frameRate,movTitle,colbarTitle)
        %%%
            % sig2Movie: channels x timepoints
            % chanMap: 2D channel map
            % timeAll: 1 x timepoints (in seconds)
            % etw: epoch time window in seconds (e.g. [-1 1] to print movie
            % between -1 to 1 seconds
            % clim: Colorbar range in uV or z-score value (e.g. [0 20])
            % frameRate: Frame rate of the movie (e.g. 120)
            % movTitle: Filename to be saved (e.g. 'S23_highGamma')
            % colbarTitle: Color axis label (e.g. 'z-score')


%timeAll = linspace(tw(1),tw(2),(tw(2)-tw(1))*fs);
        selectedChannels = sort(chanMap(~isnan(chanMap)))';
        timeSelectInd = timeAll>=etw(1)&timeAll<=etw(2);
        timeSelect = timeAll(timeSelectInd);
        figure;
        plot(timeSelect,sig2Movie(:,timeSelectInd),'color',[0 0 0] +0.75);
        hold on;
        plot(timeSelect,mean(sig2Movie(:,timeSelectInd),1),'color',[0 0 0]);
        ylim(clim);
        xlabel('Time (s)');
        ylabel(colbarTitle);
        title(movTitle);
        saveas(gcf,[movTitle '_timeSeries.png']);
        waveChannel = nan(size(chanMap,1),size(chanMap,2),length(timeAll(timeSelectInd)));
        for c = 1 : length(selectedChannels)
            [cIndR, cIndC] = find(ismember(chanMap,selectedChannels(c)));
            waveChannel(cIndR,cIndC,:)=sig2Movie(c,timeSelectInd);
        end
        figure;
        for iTime=1:size(waveChannel,3)
          %  surfc(X,Y,sq(spec_chansBHG(:,:,iT)),'FaceAlpha',0.5);
            b = imagesc(sq(waveChannel(:,:,iTime))); 
            cb = colorbar;
            ylabel(cb,colbarTitle)
            %truesize(gcf,[1000 500]);
            set(b,'AlphaData',~isnan(sq(waveChannel(:,:,iTime))));
            caxis([clim(1) clim(2)])          
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])        
            axis equal
            axis tight
            set(gca,'FontSize',20);
            colormap(parula(4096))
            
         title([num2str(round(timeSelect(iTime),2)) ' s'])
           M(iTime)=getframe(gcf);
        end
        cmap=colormap('jet');
        close
        vname = strcat(movTitle,'.avi');
        vidObj=VideoWriter(vname, 'Motion JPEG AVI');
        vidObj.Quality = 100;    
        vidObj.FrameRate = frameRate;
        open(vidObj);        
        writeVideo(vidObj,M);
         close(vidObj);
end