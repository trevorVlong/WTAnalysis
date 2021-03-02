% a design study where I compare the old and the new
close all

% load paths and dependencies
load('options.mat');

% add paths to current workspace
fdnames = fieldnames(paths);
for fdn = 1:numel(fdnames)
    addpath(paths.(fdnames{fdn}));
end

datapath = strcat(paths.dropboxpath,'/','spring21/OldData/');
addpath(datapath);

%% 
plotops.interp = 1;
plotops.dcjmax = 6.5;
plotops.colormap = turbo(1000);
plotops.alpharange = [0,25];
plotops.clrange    = [-2,10];

%% --------------------------------------------------
% pick files of interest

Re = 'H';
df = [20,40,50,55];
geom = {"solid","slotted"};
plotops.names = {"single-element","slotted"};

n = numel(geom);
m = numel(df);

files = cell(n,m);
for geomiN = 1:n
    geomi = geom{geomiN};
    for dfiM = 1:m
        dfi = df(dfiM);
        
        filename = sprintf("%s_%s_%2.0f.mat",geomi,Re,dfi);
        file_i = strcat(datapath,filename);
        
        
        try
            isfile(file_i);
            files{geomiN,dfiM} = filename;
        catch e 
            warning('file %s does not exist',filename);
            
        end
    end
end


% create grid
motormaps = {'F40II-5in','F40II-4in'};


fprintf('plot grid is %1.0f x %1.0f\n\n',n,m);


clplots = figure('Visible','off');
cxplots = figure('Visible','off');
cmplots = figure('Visible','off');
polarplots = figure('Visible','off');
%--------------------------------------------------
% new prop models
pos = 1;
for filen = 1:n
    for filem = 1:m
        
        file = files{filen,filem};
        plotops.name = plotops.names{filen};
        % plot clalpha
        set(0,'CurrentFigure',clplots)
        subplot(n,m,pos);
        plotcla(file,clplots,plotops);
        
        % plot cxalpha
        set(0,'CurrentFigure',cxplots)
        subplot(n,m,pos);
        plotcxa(file,cxplots,plotops);
        
        % plot polar
        set(0,'CurrentFigure',polarplots)
        subplot(n,m,pos);
        plotpolar(file,polarplots,plotops);
        
        % plot cm-alpha
        set(0,'CurrentFigure',cmplots)
        subplot(n,m,pos);
        plotcma(file,cmplots,plotops);        
        pos = pos+1;
    end
end

% clplots.Visible = 'on';
% cxplots.Visible = 'on';
% polarplots.Visible = 'on';
cmplots.Visible = 'on';





