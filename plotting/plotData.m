function plotData(flapN,flapT,Re,plotops)
    
    
    plotops.linewidth = 1.5;
    % default plotting options if no plotvars
    if nargin < 4
        plotops.Msize = 30;
        plotops.xlabel = "$c_x$";
        plotops.ylabel = "$c_l$";
        plotops.suptitle = "";
        plotops.fsize  = 14;
    end
    
    % hard set dcj limits
    
    % set this to dynamically change with each plot
    
    
%==========================================================================
% loop to build plot table

    for NF = 1:length(flapN)
        for FT = 1:length(flapT)
            % set index and load file
            idx = length(flapN)*(FT-1) + NF;
            file = sprintf("%s_%s_%2.0f",flapT(FT),Re,flapN(NF));
            load(file);
            %generate subplot and scatter based on plot type in plotvars
            subplot(length(flapT),length(flapN),idx);
            plotops.clim   = [0,max(dcj(:))];
            map = parula(6);
            
            if plotops.type == "polar"
                h(idx) = scatter(cx2D(:),cl2D(:),plotops.Msize,dcj(:),'filled');
                hold on
                alfamax = max(alfa,[],'all');
                alfamin = min(alfa,[],'all');
                alfa_lin = linspace(alfamin,alfamax,40);
                
                dcjvec = linspace(0,max(dcj(:)),6);
                for ndcj = 1:6
                    dcj = dcjvec(ndcj);
                    cl = Vcl(alfa_lin,dcj*ones(size(alfa_lin)));
                    cx = Vcx(alfa_lin,dcj*ones(size(alfa_lin)));
                    plot(cx,cl,'-',...
                        'Color',map(ndcj,:) ,...
                        'MarkerEdgeColor',map(ndcj,:) ,...
                        'LineWidth',plotops.linewidth);
                end
                
            elseif plotops.type == "cla"
                h(idx) = scatter(alfa(:),cl2D(:),plotops.Msize,dcj(:),'filled');
                                hold on
                alfamax = max(alfa,[],'all');
                alfamin = min(alfa,[],'all');
                alfa_lin = linspace(alfamin,alfamax,40);
                dcjvec = linspace(0,max(dcj(:)),6);
                for ndcj = 1:6
                    dcj = dcjvec(ndcj);
                    cl = Vcl(alfa_lin,dcj*ones(size(alfa_lin)));
                    plot(alfa_lin,cl,'-',...
                        'Color',map(ndcj,:) ,...
                        'MarkerEdgeColor',map(ndcj,:) ,...
                        'LineWidth',plotops.linewidth);
                end
                
            elseif plotops.type == "cxa"
                h(idx) = scatter(alfa(:),cx2D(:),plotops.Msize,dcj(:),'filled');
                                hold on
                alfamax = max(alfa,[],'all');
                alfamin = min(alfa,[],'all');
                alfa_lin = linspace(alfamin,alfamax,40);
                dcjvec = linspace(0,max(dcj(:)),6);
                for ndcj = 1:6
                    dcj = dcjvec(ndcj);
                    cx = Vcx(alfa_lin,dcj*ones(size(alfa_lin)));
                    plot(alfa_lin,cx,'-',...
                        'Color',map(ndcj,:) ,...
                        'MarkerEdgeColor',map(ndcj,:) ,...
                        'LineWidth',plotops.linewidth);
                end
                
            elseif plotops.type == "cma"
                h(idx) = scatter(alfa(:),cm2D(:),plotops.Msize,dcj(:),'filled');
                                hold on
                alfamax = max(alfa,[],'all');
                alfamin = min(alfa,[],'all');
                alfa_lin = linspace(alfamin,alfamax,40);
                dcjvec = linspace(0,max(dcj(:)),6);
                for ndcj = 1:6
                    dcj = dcjvec(ndcj);
                    cm = Vcm(alfa_lin,dcj*ones(size(alfa_lin)));
                    plot(alfa_lin,cm,'-',...
                        'Color',map(ndcj,:) ,...
                        'MarkerEdgeColor',map(ndcj,:) ,...
                        'LineWidth',plotops.linewidth);
                end
            end
            
            
            %-------------------------------------------------
            %labels and limits for each subplot
            xlim(plotops.xlim);
            xlabel(plotops.xlabel,'interpreter','latex','FontSize',plotops.fsize);
            
            ylim(plotops.ylim);
            ylabel(plotops.ylabel,'interpreter','latex','FontSize',plotops.fsize);
            
            title(sprintf("%s $\\delta_f$ = %2.0f$^{\\circ}$",plotops.names(FT),flapN(NF)),'interpreter','latex')
            set(gca,'FontName','CMU Serif')
            grid on

        end
    end
    
    %======================================================================
    % polishing plot
    
    % colorbar
    colormap('parula')
    caxis(plotops.clim)
    cb = colorbar('Position',[.935,.1,.02,.8]);
    cb.Label.String = "\Delta CJ";
    cb.Label.FontSize = plotops.fsize;
    cb.Label.Rotation = 0;
    
    
    % chart title and macros
    sgtitle(plotops.suptitle)
    set(gca,'FontName','CMU Serif')
    plotops.fh.WindowState = 'maximized';
    pause(2);
    
    %optional save
    if plotops.save == "save"
        saveas(gcf,plotops.plotname, 'epsc');
    end
end