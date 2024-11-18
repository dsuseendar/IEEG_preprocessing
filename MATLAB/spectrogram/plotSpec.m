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
        args.cval = [0.7 1.4];
    end

    tspec = linspace(tw(1), tw(2), size(spec, 1)); % Time vector for spectrograms
    if(args.isdb)
        spec2plot = 20.* log10(squeeze(spec));
    else
        spec2plot = squeeze(spec);
    end
    imagesc(tspec, [], spec2plot');
    xlim(args.etw);
    set(gca,'YTick',1:args.f_interal:length(freq));
    set(gca,'YTickLabels',round(freq(1:args.f_interal:end)));
    set(gca,'FontSize',args.font_size);
    title(args.title);
    if(args.db)
        caxis(20.*log10(cval))
    else
        caxis(cval)
    end
end