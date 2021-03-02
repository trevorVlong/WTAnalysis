close all


alfarange = 0:3:15;
dcjrange  = 0:0.5:8;
dfrange   = 0:10:50;

load('options.mat','plotops','paths');
plotops.colormap = turbo(1000);
plotops.dcjmax = 8;
plotops.interp = 1;
plotops.analytic = 0;


% setup figures
clafig = figure('Visible','off');
cxafig = figure('Visible','off');
cmafig = figure('Visible','off');
polarfig = figure('Visible','off');

%% Plot Spence Model





hold on;
df = 20;

if plotops.analytic == 1
    for dcjN = 1:4:length(dcjrange)
        dcjval = dcjrange(dcjN);
        cl = [];
        for alfaval = alfarange
            valtable = SpenceModel(dcjval,alfaval,df,3);
            cl = [cl, valtable.cl];
        end

        figure(clafig);
        hold on;
        ax = plot(alfarange,cl,...
             'LineWidth',1.5,...
             'Color',plotops.colormap(dcjN,:));     
    end
    

    
end
    

    
%% plot results
geom = 'solid_mod';
Re   = 'L';
flapd= df;

filename = sprintf('%s_%s_%02.0f.mat',geom,Re,flapd);

plotcla(filename,clafig,plotops);
plotcxa(filename,cxafig,plotops);
plotpolar(filename,polarfig,plotops);


% -----------------------------------------------
% turn figures back on

clafig.Visible = 'on';
cxafig.Visible = 'on';
% cmafig.Visible = 'on';
polarfig.Visible = 'on';
    

