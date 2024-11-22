function plotSpec(spec,tw,freq,args)
    arguments
        spec
        tw
        freq
        args.isdb = true;
        args.f_interal = 5;
        args.font_size = 8;
        args.etw = tw;
        args.title = '';
        args.cval = [-2 2];
        
    end

    tspec = linspace(tw(1), tw(2), size(spec, 1)); % Time vector for spectrograms
    
    imagesc(tspec, [], spec2plot');
    
    set(gca,'YTick',1:args.f_interal:length(freq));
    set(gca,'YTickLabels',round(freq(1:args.f_interal:end)));
    set(gca,'FontSize',args.font_size);
    set(gca, 'YDir', 'normal');
    xlim(args.etw);
   
    title(args.title,'Interpreter','none');
    caxis(args.cval)
    
end