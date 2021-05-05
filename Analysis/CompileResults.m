% plots results from a single set of runs on tiled plots

clear all
setup_analysis()

calnum  = [2];
% tarenum = [124:132];
% tarenum = [133:140]
%tarenum = [141:149]; % lowered mounts
tarenum = 115:123; %nominal mounts
%tarenum = 133:140; % slotted flaps
% tarenum = 124:132; %small prop
%tarenum = [113,112,108,104,105,106,109,110,111]; % high mount
% tarenum = 141:149; % low mount
% tarenum = 150:158; % 20 deg mount
%tarenum = 159:167; % 0  deg mount
% tarenum = [113,112,108,105,106,107,109,110,111]; high mount
dfvals = [0,20,30,35,40,45,50,55,60];
% tarenum = 131;
runnum = 1:11;
% runnum = 6;

config = 'BaseConfig';
plotting = 1;
saverun = 1;


datadir = "/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/spring21/PerformanceRuns/";
rundir  = strcat(datadir,"runs/runs");
caldir  = strcat(datadir,"calibrations");
taredir = strcat(datadir,"tares");
addpath(datadir,rundir,taredir,caldir);

load('analysis_options.mat','statops','plotops');

if statops.rmverbose == 1
   FID = fopen('log.txt','w'); 
end
% set collection matrices
perftable = struct; 

alfaoffset = 3.7;


calfiles = {};
tarefiles = {};
runfiles = {};
runtable = {};
cl = [];
cx = [];
cm = [];
std_cl = [];
std_cx = [];
std_cm = [];
min_cl = [];
min_cx = [];
min_cm = [];
max_cl = [];
max_cx = [];
max_cm = [];

errorfiles = [];

dcj = [];
maxdcj = 0;
df  = [];
alfa = [];
cidx = 1;
for cnum = calnum
    calfile = sprintf("cal%02.0f.mat",cnum);
    
    if isfile(strcat(caldir,'/',calfile))
        calfiles{cnum} = calfile;
        cidx = cidx+1;
    end
    
    tidx = 0;
    for tnum = tarenum
        % each flap position
        
        tarefile = sprintf("tare%03.0f.mat",tnum);
        
        if isfile(strcat(taredir,'/',tarefile))
            tidx = tidx+1;
            tarefiles{tidx} = tarefile;
        end
        
        
        for rnum = runnum
            % each angle of attack
            n = 1;
            
            
            for runpos = 1:20
                % each power setting
                runfile = sprintf("run_%02.0f_%03.0f_%03.0f_%02.0f.mat",cnum,tnum,rnum,runpos);
                if isfile (strcat(rundir,'/',runfile))
                    try
                        runfiles{rnum,n,tidx} = runfile;
                        
                        runtable{rnum,n,tidx} = nonDim(runfile,calfile,tarefile,FID);
                    catch e 
                        wng = sprintf("Problem analyzing %s, check file directly",runfile);
                        warning(wng);
                        errorfiles = [errorfiles;runfile];
                        break
                    end
                    
                    %------------------------------------------------------
                    % values of interest to save in data structures
                    
                    %lift
                    cl(n,rnum,tidx) = runtable{rnum,n,tidx}.cl_average;
                    std_cl(n,rnum,tidx) = runtable{rnum,n,tidx}.cl_std;
                    min_cl(n,rnum,tidx) = runtable{rnum,n,tidx}.min_cl;
                    max_cl(n,rnum,tidx) = runtable{rnum,n,tidx}.max_cl;
                    

                    % x-force
                    Xforce(n,rnum,tidx) = runtable{rnum,n,tidx}.Fx;
                    cx(n,rnum,tidx) = runtable{rnum,n,tidx}.cx_average;
                    std_cx(n,rnum,tidx) = runtable{rnum,n,tidx}.cx_average;
                    min_cx(n,rnum,tidx) = runtable{rnum,n,tidx}.cx_average;
                    max_cx(n,rnum,tidx) = runtable{rnum,n,tidx}.cx_average;
                    
                    % drag
                    Drag(n,rnum,tidx) = runtable{rnum,n,tidx}.Fd;
                    cd(n,rnum,tidx) = runtable{rnum,n,tidx}.cd_average;
                    
                    
                    % moment
                    cm(n,rnum,tidx) = runtable{rnum,n,tidx}.cm_average;
                    std_cm(n,rnum,tidx) = runtable{rnum,n,tidx}.cm_average;
                    min_cm(n,rnum,tidx) = runtable{rnum,n,tidx}.cm_average;
                    max_cm(n,rnum,tidx) = runtable{rnum,n,tidx}.cm_average;
                    
                    %------------------------------------------------------
                    % jet performance
                    dcj(n,rnum,tidx) = runtable{rnum,n,tidx}.dCJ;
                    std_dcj(n,rnum,tidx) = runtable{rnum,n,tidx}.dCJstd;
                    min_dcj(n,rnum,tidx) = runtable{rnum,n,tidx}.mindCJ;
                    max_dcj(n,rnum,tidx) = runtable{rnum,n,tidx}.maxdCJ;
                    
                    if dcj(n,rnum,tidx) > maxdcj
                        maxdcj = dcj(n,rnum,tidx);
                    end
                    
                    % CT,RPM and std for each motor for diagnostic
                    rpm1(n,rnum,tidx) = runtable{rnum,n,tidx}.rpm1;
                    rpm2(n,rnum,tidx) = runtable{rnum,n,tidx}.rpm2; 
                    rpm3(n,rnum,tidx) = runtable{rnum,n,tidx}.rpm3;
                    rpm4(n,rnum,tidx) = runtable{rnum,n,tidx}.rpm4;
                    
                    CT(n,rnum,tidx) = mean(runtable{rnum,n,tidx}.CT);
                    perftable.Thrust(n,rnum,tidx) = sum(runtable{rnum,n,tidx}.Thrust);
                 
                    std_rpm1(n,rnum,tidx) = runtable{rnum,n,tidx}.RPMstd(1);
                    std_rpm2(n,rnum,tidx) = runtable{rnum,n,tidx}.RPMstd(2);
                    std_rpm3(n,rnum,tidx) = runtable{rnum,n,tidx}.RPMstd(3);
                    std_rpm4(n,rnum,tidx) = runtable{rnum,n,tidx}.RPMstd(4);
                    
                    min_rpm1(n,rnum,tidx) = runtable{rnum,n,tidx}.minRPM(1);
                    min_rpm2(n,rnum,tidx) = runtable{rnum,n,tidx}.minRPM(2);
                    min_rpm3(n,rnum,tidx) = runtable{rnum,n,tidx}.minRPM(3);
                    min_rpm4(n,rnum,tidx) = runtable{rnum,n,tidx}.minRPM(4);
                    
                    max_rpm1(n,rnum,tidx) = runtable{rnum,n,tidx}.maxRPM(1);
                    max_rpm2(n,rnum,tidx) = runtable{rnum,n,tidx}.maxRPM(2);
                    max_rpm3(n,rnum,tidx) = runtable{rnum,n,tidx}.maxRPM(3);
                    max_rpm4(n,rnum,tidx) = runtable{rnum,n,tidx}.maxRPM(4);
                    
                    
                    % geom and other
                    df(n,rnum,tidx)  = runtable{rnum,n,tidx}.flap;
                    alfa(n,rnum,tidx) = mean(runtable{rnum,n,tidx}.alfa) + alfaoffset;
                    Vinf(n,rnum,tidx) = mean(runtable{rnum,n,tidx}.Vinf);
                    rho(n,rnum,tidx) = mean(runtable{rnum,n,tidx}.rho);
                    
                    
                    n = n+1;
                end

                clear runfile
            % for each power
            end
            % for each AoA
        end
        alf = alfa(:,:,tidx);
        dc  = dcj(:,:,tidx);
        clv = cl(:,:,tidx);
        cxv = cx(:,:,tidx);
        cmv = cm(:,:,tidx);
        perftable.Vcla{tidx} = scatteredInterpolant(alf(:),dc(:),clv(:),'linear','none');
        perftable.Vcxa{tidx} = scatteredInterpolant(alf(:),dc(:),cxv(:),'linear','none');
        perftable.Vcma{tidx} = scatteredInterpolant(alf(:),dc(:),cmv(:),'linear','none');
        % for each flap setting
    end
end

fclose(FID);
%% Data handling and analysis
% assign all data to table


% aero perforamnce
perftable.calfiles = calfiles;
perftable.tarefiles = tarefiles;
perftable.runfiles = runfiles;
perftable.runtable = runtable;
perftable.cl = cl;
perftable.cx = cx;
perftable.cd = cd;
perftable.cm = cm;

perftable.Xforce = Xforce;
perftable.Drag = Drag;


perftable.std_cl = std_cl;
perftable.std_cx = std_cx;
perftable.std_cm = std_cm;
perftable.min_cl = min_cl;
perftable.min_cx = min_cx;
perftable.min_cm = min_cm;
perftable.max_cl = max_cl;
perftable.max_cx = max_cx;
perftable.max_cm = max_cm;


% motors and thrust
perftable.CT = CT;

perftable.rpm1 = rpm1;
perftable.rpm2 = rpm2;
perftable.rpm3 = rpm3;
perftable.rpm4 = rpm4;

perftable.min_rpm1 = min_rpm1;
perftable.min_rpm2 = min_rpm2;
perftable.min_rpm3 = min_rpm3;
perftable.min_rpm4 = min_rpm4;

perftable.max_rpm1 = max_rpm1;
perftable.max_rpm2 = max_rpm2;
perftable.max_rpm3 = max_rpm3;
perftable.max_rpm4 = max_rpm4;

perftable.errorfiles = errorfiles;

perftable.dcj = dcj;
perftable.dcjcolor = [];
perftable.df  = df;
perftable.alfa = alfa;
perftable.rho = rho;
perftable.Vinf = Vinf;



%% setup for plotting

% remove nans
perftable.cl(perftable.cl == 0) = 999;
perftable.cx(perftable.cx == 0) = 999;
perftable.cm(perftable.cm == 0) = 999;
perftable.alfa(perftable.alfa == 0) = 999;

perftable.dcj(perftable.dcj == 0) = 999;
perftable.dcj(1,:,:) = 0;


% get dcj color index / rgb triplet (use k steps)
% add to plotting options 
k = 1000;
cmap = turbo(k);
dcjcolor_index = floor((perftable.dcj/(maxdcj + 1))*1000)+1;

% for idx = 1:length(dcj(:))
%     try
%         perftable.dcjcolor(idx,:) = cmap(dcjcolor_index(idx),:);
%     catch
%         perftable.dcjcolor(idx,:) = [NaN, NaN, NaN];
%     end
% end



%%
if plotting == 1
    close all
    % interpolant
    clalfafig = figure('Visible','off');
    cxalfafig = figure('Visible','off');
    cdalfafig = figure('Visible','off');
    cmalfafig = figure('Visible','off');
    polarfig = figure('Visible','off');


    for flapnum = 1:length(perftable.cl(1,1,:))

        clflap = perftable.cl(:,:,flapnum);
        clflap(clflap == 0) = nan;

        cxflap = perftable.cx(:,:,flapnum);
        cxflap(cxflap == 0) = nan;
        
        cdflap = perftable.cd(:,:,flapnum);
        cd(cd == 0) = nan;

        cmflap = perftable.cm(:,:,flapnum);
        cmflap(cmflap == 0) = nan;

        alfaflap = perftable.alfa(:,:,flapnum);
        alfaflap(alfaflap == 0) = nan;

        dcjflap = perftable.dcj(:,:,flapnum);
        dcjflap(dcjflap == 999) = nan;

        % cl-alpha
        set(0,'CurrentFigure',clalfafig);
        subplot(3,3,flapnum)
    %     [out] = plotcl-alpha()
        scatter(alfaflap(:),clflap(:),40,dcjflap(:),'filled');
        ylim([-1,10]);
        xlim([-15,30]);
        title(sprintf('$\\delta_f$ = %2.0f',dfvals(flapnum)),'Interpreter','latex');
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('$c_\ell$','Interpreter','latex')
        grid on
        colorbar
        ax = gca;
        set(ax,'FontSize',13)
        set(ax,'FontName','avenir')    

        % cm- alpha

        set(0,'CurrentFigure',cmalfafig);
        subplot(3,3,flapnum)
        scatter(alfaflap(:),cmflap(:),40,dcjflap(:),'filled');
        ylim([-1,1]);
        xlim([-15,30]);
        title(sprintf('$\\delta_f$ = %2.0f',dfvals(flapnum)),'Interpreter','latex');
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('$c_m$','Interpreter','latex')
        grid on
        colorbar
        ax = gca;
        set(ax,'FontSize',13)
        set(ax,'FontName','avenir')    


        % cx- alpha

        set(0,'CurrentFigure',cxalfafig);
        subplot(3,3,flapnum)
        scatter(alfaflap(:),cxflap(:),40,dcjflap(:),'filled');
        ylim([-4,1]);
        xlim([-15,30]);
        title(sprintf('$\\delta_f$ = %2.0f',dfvals(flapnum)),'Interpreter','latex');
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('$c_x$','Interpreter','latex')
        grid on
        colorbar
        ax = gca;
        set(ax,'FontSize',13)
        set(ax,'FontName','avenir')

        % cl- cx

        set(0,'CurrentFigure',polarfig);
        subplot(3,3,flapnum)
        scatter(cxflap(:),clflap(:),40,dcjflap(:),'filled');
        ylim([-1,10]);
        xlim([-4,1]);
        title(sprintf('$\\delta_f$ = %2.0f',dfvals(flapnum)),'Interpreter','latex');
        xlabel('$c_x$','Interpreter','latex')
        ylabel('$c_\ell$','Interpreter','latex')
        grid on
        colorbar
        ax = gca;
        set(ax,'FontSize',13)
        set(ax,'FontName','avenir')
        % other?
        
        
        
        set(0,'CurrentFigure',cdalfafig);
        subplot(3,3,flapnum)
        scatter(alfaflap(:),cdflap(:),40,dcjflap(:),'filled');
        ylim([-1,5]);
        xlim([-10,30]);
        title(sprintf('$\\delta_f$ = %2.0f',dfvals(flapnum)),'Interpreter','latex');
        xlabel('$\alpha$','Interpreter','latex')
        ylabel('$c_d$','Interpreter','latex')
        grid on
        colorbar
        ax = gca;
        set(ax,'FontSize',13)
        set(ax,'FontName','avenir')
        % other?
        
        
    end
    
% turn figures back on
    clalfafig.WindowState = "maximized";
    sgtitle(clalfafig,sprintf('c_l - \\alpha %s',config))
    clalfafig.Visible = 'on';

    cmalfafig.WindowState = "maximized";
    sgtitle(cmalfafig,sprintf('c_m - \\alpha %s',config))
    cmalfafig.Visible = 'on';

    cxalfafig.WindowState = "maximized";
    sgtitle(cxalfafig,sprintf('c_x - \\alpha %s',config))
    cxalfafig.Visible = 'on';

    polarfig.WindowState = "maximized";
    sgtitle(polarfig,sprintf('c_l-c_x %s',config))
    polarfig.Visible = 'on';
    
    cdalfafig.WindowState = "maximized";
    sgtitle(cdalfafig,sprintf('c_d - \\alpha %s',config))
    cdalfafig.Visible = 'on';
    % n-dimensional fit

    % other???
end

%% save data structures
if saverun == 1
    
    % data
    save(sprintf('%s.mat',config),'perftable');
    
    % figures
    if plotting == 1
        saveas(clalfafig,sprintf('cl_alfa_%s',config),plotops.savetype);
        saveas(cxalfafig,sprintf('cx_alfa_%s',config),plotops.savetype);
        saveas(cmalfafig,sprintf('cm_alfa_%s',config),plotops.savetype);
        saveas(cdalfafig,sprintf('cd_alfa_%s',config),plotops.savetype);
        saveas(polarfig,sprintf('polar_%s',config),plotops.savetype)
    end
end
