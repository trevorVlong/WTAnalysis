% a design study where I compare the old and the new
close all

% load paths and dependencies
load('options.mat');

% add paths to current workspace
fdnames = fieldnames(paths);
for fdn = 1:numel(fdnames)
    addpath(paths.(fdnames{fdn}));
end

datapath = strcat(paths.dropboxpath,'/','sum_fall_20/');
addpath(datapath);

%% 
plotops.interp = 1;
plotops.dcjmax = 6.5;
plotops.colormap = turbo(1000);
plotops.alpharange = [0,25];
plotops.clrange    = [-2,10];

%% --------------------------------------------------
% pick files of interest
runfile = 'run_10_303_016_25.mat';
calfile = 'cal10.mat';
tarefile = 'tare303.mat';

load(runfile);


wallfig = figure('Visible','off');
wd = mean(walldata([1:11,13:29,12],:),2);
wdc = wd(1:28)/wd(end)-1;
plot(wdc,'b',...
     'LineWidth',2)
grid on
set(gca,'FontName',plotops.font);
set(gca,'FontSize',10);
xlabel('pitot position number','Interpreter','latex',...
        'FontSize',12);
ylabel('$\frac{p_{wall,static}}{p_\infty}$','Interpreter','latex',...
        'Rotation',0,...
        'Position',[-3,.95],...
        'FontSize',20);
set(wallfig,'Visible','on');
%% 
[corrected_data,uncorrected_data] = nonDimVecs(runfile,calfile,tarefile);

clu = uncorrected_data.cl_average;
cxu = uncorrected_data.cx_average;
Vinf = uncorrected_data.Vinf;

clv = corrected_data.clv;
cxv = corrected_data.cxv;
Veff= corrected_data.Veffv;
alfa= corrected_data.alfa;

clvec = [clu clv];
cxvec = [cxu cxv];
Veffvec = [Vinf Veff];


% other information
dcj = corrected_data.dCJ;
df  = corrected_data.flap;
convplot = figure('Visible','off');
hold on;

plot(clvec,'b',...
    'LineWidth',2)
plot(clvec(1)*ones(size(clvec)),'r--');
plot(clvec(end)*ones(size(clvec)),'r--');


set(gca,'FontName',plotops.font);
set(gca,'FontSize',10);
xlabel('Iteration Number','Interpreter','latex',...
        'FontSize',12);
ylabel('$c_\ell$ corrected','Interpreter','latex')
title(sprintf('$\\delta_f = %2.0f$, $\\alpha = %2.1f$ , $\\Delta{c_J} = %2.2f$',df,alfa,dcj),'Interpreter','latex');

set(convplot,'Visible','on')