function cplot = wakeContour(uVdata,runtable,plotops)
% plots contour plot of wake for 3x2 wind tunnel runs. Does NOT doe labels
% and titles, those must be done outside of the loop

    try
        [a,cplot] = contourf(linspace(0,runtable.Nspan,runtable.Nspan)/24,runtable.zpos/9,...
        round(uVdata,2),...
                'LineWidth',1,...
                'LevelList',plotops.VJ_rng); 
        caxis([plotops.VJ_rng(1),plotops.VJ_rng(end)]);
        
        % xlabel and ylabel
        xlbl = xlabel("span position X/b",'interpreter','latex');
        ylbl = ylabel("vertical position Z/c",'interpreter','latex');

        xlbl.FontSize = 13;
        ylbl.FontSize = 13;
    catch
        warning('something went wrong plotting')
        cplot = contourf(NaN,NaN,NaN);
    end

end