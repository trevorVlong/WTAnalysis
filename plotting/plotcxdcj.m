function [] = plotcxdcj(file,fig,plotops)
%PLOTCLA plots cl vs alpha with dcj sidebar
% datatable takes in a table 

% ===========================================
% load runfile

    load(file);
    

% ----------------------------------- 
    set(0,'CurrentFigure',fig)
    hold on
    colormap(plotops.colormap);
    ax = scatter(dcj(:),cx2D(:),20,'filled');
    
    

    
    if plotops.interp == 1
        alfarange = [min(alfa(:)), max(alfa(:))];
        dcjrange  = [0, max(dcj(:))];
        
        alfavec = linspace(alfarange(1),alfarange(2),20);
        dcjvec  = linspace(0,dcjrange(2),5);
        
        for ndcj = 1:5
            dcjcolorindex = round(dcjvec(ndcj)/plotops.dcjmax*1000);
            
            clinterp = Vcl(alfavec,dcjvec(ndcj)*ones(size(alfavec)));
            
            plot(alfavec,clinterp,...
                 'Color', plotops.colormap(dcjcolorindex+1,:),...
                 'LineWidth',1.5 );
            hold on
            
        end
        
        
    end
    
    % colorbar etc
    cb = colorbar;
    caxis([0,plotops.dcjmax]);
    cb.Label.Interpreter= 'latex';
    cb.Label.String = "$\Delta{c_J}$";
    cb.Label.Rotation = 0;
    cb.Label.Position = [3,4.75,0];
    
    % fontsize
    set(gca,'FontName',plotops.font);
    set(gca,'FontSize',10);
    grid on   
    
    xlim(plotops.dcjrange);
    ylim(plotops.clrange);
    % labels, 
    xlabel('$\alpha$','Interpreter','latex','FontSize',plotops.labelFS);
    ylabel('$c_\ell$','Interpreter','latex','FontSize',plotops.labelFS);
    cb.Label.FontSize = plotops.labelFS;
    
    ttl = sprintf('%s $\\delta_f$ = %2.0f',plotops.name, flap_ang);
    title(ttl,'Interpreter','latex');

end

