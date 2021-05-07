% overlays performance plots of two different geometires by specifying the
% file names

compfiles = ["BaseConfig","HighMount","LowMount"];
leg       = ["Base Configuration","Raised Mount","Lowered Mount"];

scatterpoints = 0;
interpolants = 1;

show_polar = 1;
show_cm_cl = 1;
show_falfa = 1;


savecomp = 1;

% flap settings to compare
compare =   [1, 1, 1, 1, 1, 1, 1, 1, 1 ];
df=         [0  20 30 35 40 45 50 55 60];
len = length(df);

dcjq = [1,2,3,4,4.5,5];
alfaq = linspace(-15,30,10);

% load any needed options for plotting
addpath('../');
load('analysis_options.mat','plotops');
ms = 50;


% figures
figs = gobjects(9,1);
polarfigs = gobjects(9,1);
fclfigs = gobjects(9,1);
plotcolumns = 1;


for idx = 1:len
    
    if compare(idx) == 1
        figs(idx) = figure('Visible','off');
        for file_idx = 1:length(compfiles)
            
            
            load(compfiles(file_idx));
            comptable = perftable;
            
            if compfiles(file_idx) == "SlottedFlap"
                idx = idx-1;
            end
            Vcl = comptable.Vcla{idx};
            Vcx = comptable.Vcxa{idx};
            Vcm = comptable.Vcma{idx};
            clcomp = comptable.cl(:,:,idx);
            clcomp(clcomp > 20) = nan;
            cxcomp = comptable.cx(:,:,idx);
            cxcomp(abs(cxcomp) > 20) = nan;
            cmcomp = comptable.cm(:,:,idx);
            cmcomp(cmcomp > 20) = nan;
            cdcomp = comptable.cd(:,:,idx);
            cdcomp(abs(cdcomp) > 20) = nan;
            alfacomp = comptable.alfa(:,:,idx);
            alfacomp(alfacomp > 40) = nan;
            dcjcomp = comptable.dcj(:,:,idx);
            dcjcomp(dcjcomp > 20) = nan;
            
           %clalpha
            h = subplot(3,plotcolumns,1);
            
            if file_idx == 1
                if scatterpoints == 1
                    scatter(alfacomp(:),clcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    outdata = plotinterpolants(h,alfaq,dcjq,Vcl,file_idx);
                end
                
                colorbar
                caxis([0,6]);
                grid on
                grid minor
                ylabel('$c_l$','Interpreter','latex')
                ylim([-1,10]);
            else
%                 scatter(alfacomp(:),clcomp(:),ms,dcjcomp(:),...
%                     'Marker',plotops.markerstyle{file_idx}); hold on;

                if scatterpoints == 1
                    scatter(alfacomp(:),clcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    outdata = plotinterpolants(h,alfaq,dcjq,Vcl,file_idx);
                end
                if file_idx == length(compfiles)
                    set(gca,'FontSize',plotops.axesFS)
                end
            end


           %cxalpha
            h = subplot(3,plotcolumns,2);
            
            if file_idx == 1
                if scatterpoints == 1
                    scatter(alfacomp(:),cxcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    outdata = plotinterpolants(h,alfaq,dcjq,Vcx,file_idx);
                end
                
                colorbar
                caxis([0,6]);
                grid on
                grid minor
                ylabel('$c_x$','Interpreter','latex')
                ylim([-3,1]);
            else
                if scatterpoints == 1
                    scatter(alfacomp(:),cxcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    outdata = plotinterpolants(h,alfaq,dcjq,Vcx,file_idx);
                end
                if file_idx == length(compfiles)
                    set(gca,'FontSize',plotops.axesFS)
                end
            end
            
            
           %cmalpha
            h = subplot(3,plotcolumns,3);
            
            if file_idx == 1
                if scatterpoints == 1
                    scatter(alfacomp(:),cmcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    outdata = plotinterpolants(h,alfaq,dcjq,Vcm,file_idx);
                end
                
            colorbar
            caxis([0,6]);
            grid on
            grid minor
            xlabel('$\alpha$','Interpreter','latex')
            ylabel('$c_m$','Interpreter','latex')
            ylim([-1,0.5]);
            else
                if scatterpoints == 1
                    scatter(alfacomp(:),cmcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    outdata = plotinterpolants(h,alfaq,dcjq,Vcm,file_idx);
                end
                
                if file_idx == length(compfiles)
                    set(gca,'FontSize',plotops.axesFS)
                end
            end
        end
    end 
end

for idx = 1:len
    
    if compare(idx) == 1
        polarfigs(idx) = figure('Visible','off');
        for file_idx = 1:length(compfiles)
            
            
            load(compfiles(file_idx));
            comptable = perftable;
            
            if compfiles(file_idx) == "SlottedFlap"
                idx = idx-1;
            end
            Vcl = comptable.Vcla{idx};
            Vcx = comptable.Vcxa{idx};
            Vcm = comptable.Vcma{idx};
            clcomp = comptable.cl(:,:,idx);
            clcomp(clcomp > 20) = nan;
            cxcomp = comptable.cx(:,:,idx);
            cxcomp(abs(cxcomp) > 20) = nan;
            cmcomp = comptable.cm(:,:,idx);
            cmcomp(cmcomp > 20) = nan;
            cdcomp = comptable.cd(:,:,idx);
            cdcomp(abs(cdcomp) > 20) = nan;
            alfacomp = comptable.alfa(:,:,idx);
            alfacomp(alfacomp > 40) = nan;
            dcjcomp = comptable.dcj(:,:,idx);
            dcjcomp(dcjcomp > 20) = nan;
 
           %polar
            h = subplot(1,1,1);
            if file_idx == 1
                if scatterpoints == 1
                    scatter(clcomp(:),cxcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    %outdata = plotinterpolants(h,alfaq,dcjq,Vcl,file_idx);
                    [outdata1,outdata2] = plotInterpolantPolars(h,...
                                                                alfaq,...
                                                                dcjq,...
                                                                Vcl,...
                                                                Vcx,...
                                                                file_idx);
                end
                
            colorbar
            caxis([0,6]);
            grid on
            grid minor
            xlabel('$c_\ell$','Interpreter','latex')
            ylabel('$c_x$','Interpreter','latex')
            xlim([-2,10]);
            ylim([-5,1])
            else
                
                if scatterpoints == 1
                    scatter(clcomp(:),cxcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    [outdata1,outdata2] = plotInterpolantPolars(h,...
                                                                alfaq,...
                                                                dcjq,...
                                                                Vcl,...
                                                                Vcx,...
                                                                file_idx);
                end
                
                
                if file_idx == length(compfiles)
                    set(gca,'FontSize',plotops.axesFS)
                end
            end
        end
    end 
end

for idx = 1:len
    
    if compare(idx) == 1
        fclfigs(idx) = figure('Visible','off');
        for file_idx = 1:length(compfiles)
            
            
            load(compfiles(file_idx));
            comptable = perftable;
            
            if compfiles(file_idx) == "SlottedFlap"
                idx = idx-1;
            end
            
            % extract values
            Vcl = comptable.Vcla{idx};
            Vcx = comptable.Vcxa{idx};
            Vcm = comptable.Vcma{idx};
            clcomp = comptable.cl(:,:,idx);
            clcomp(clcomp > 20) = nan;
            cxcomp = comptable.cx(:,:,idx);
            cxcomp(abs(cxcomp) > 20) = nan;
            cmcomp = comptable.cm(:,:,idx);
            cmcomp(cmcomp > 20) = nan;
            cdcomp = comptable.cd(:,:,idx);
            cdcomp(abs(cdcomp) > 20) = nan;
            alfacomp = comptable.alfa(:,:,idx);
            alfacomp(alfacomp > 40) = nan;
            dcjcomp = comptable.dcj(:,:,idx);
            dcjcomp(dcjcomp > 20) = nan;
 
           %polar
            h = subplot(1,1,1);
            if file_idx == 1
                if scatterpoints == 1
                    scatter(clcomp(:),cmcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    %outdata = plotinterpolants(h,alfaq,dcjq,Vcl,file_idx);
                    [outdata1,outdata2] = plotInterpolantPolars(h,...
                                                                alfaq,...
                                                                dcjq,...
                                                                Vcl,...
                                                                Vcm,...
                                                                file_idx);
                end
                
            colorbar
            caxis([0,6]);
            grid on
            grid minor
            xlabel('$c_\ell$','Interpreter','latex')
            ylabel('$c_m$','Interpreter','latex')
            ylim([-3,3]);
            xlim([-1,10])
            else
                if scatterpoints == 1
                    scatter(clcomp(:),cmcomp(:),ms,dcjcomp(:),'filled',...
                        'Marker',plotops.markerstyle{file_idx}); hold on;
                end
                
                if interpolants == 1
                    %outdata = plotinterpolants(h,alfaq,dcjq,Vcl,file_idx);
                    [outdata1,outdata2] = plotInterpolantPolars(h,...
                                                                alfaq,...
                                                                dcjq,...
                                                                Vcl,...
                                                                Vcm,...
                                                                file_idx);
                end
                
                if file_idx == length(compfiles)
                    set(gca,'FontSize',plotops.axesFS)
                end
            end
        end
    end 
end


for idx = 1:length(figs)
    if compare(idx) == 1
        
        if show_falfa == 1
            cfig = figs(idx);
            set(cfig,'Visible','on')
            set(cfig,'WindowState','maximized');
            % title and legend addition
            ttl = sprintf("Performance %2.0f flaps",df(idx));
            sgtitle(ttl)
            subplot(3,plotcolumns,2);
            if interpolants == 1
                for i = 1:length(compfiles)
                    plot(nan,nan,plotops.linestyle{i})
                end
            end
            legend(leg,'Location',"best")
            if savecomp == 1
                saveas(cfig,ttl,'fig');
            end            
           
        end
        
        if show_polar == 1
            cfig = polarfigs(idx);
            set(cfig,'Visible','on')
            set(cfig,'WindowState','maximized');
            % title and legend addition
            ttl = sprintf("cl-cx %2.0f flaps",df(idx));
            sgtitle(ttl);
            if interpolants == 1
                for i = 1:length(compfiles)
                    plot(nan,nan,plotops.linestyle{i})
                end
            end
            
            legend(leg,'Location',"best")
            if savecomp == 1
                saveas(cfig,ttl,'fig');
            end
            
        end
        
        if show_cm_cl == 1
            cfig = fclfigs(idx);
            set(cfig,'Visible','on')
            set(cfig,'WindowState','maximized');
            % title and legend addition
            ttl = sprintf("cl-cm %2.0f flaps",df(idx));
            sgtitle(ttl)
            if interpolants == 1
                for i = 1:length(compfiles)
                    plot(nan,nan,plotops.linestyle{i})
                end
            end
            
            legend(leg,'Location',"best")
            
            if savecomp == 1
                saveas(cfig,ttl,'fig');
            end
            
        end
        
        

    end
end
