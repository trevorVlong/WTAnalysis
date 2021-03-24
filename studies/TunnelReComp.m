% study to look at the Re of the tunnel and it's affect on data
clear 

% =========================================
% geometry etc for re 
l = convlength(9,'in','m'); % chord length of tunnel model
mu = 1.81e-5; % dynamic viscosity of air


% =========================================
calnum  = [1];
tarenum = [2,3];
runnum = 1:1:11;

homedir = "/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/spring21/PerformanceRuns/";
rundir  = strcat(homedir,"runs");
caldir  = strcat(homedir,"calibrations");
taredir = strcat(homedir,"tares");

addpath(homedir,rundir,taredir,caldir);

calfiles = {};
tarefiles = {};
runfiles = {};
runtables = {};

% check to see if files exist and add them to search arrays and get table
for cnum = calnum
    calfile = sprintf("cal%02.0f.mat",cnum);
    
    if isfile(strcat(caldir,'/',calfile))
        calfiles{cnum} = calfile;
    end
    
    for tnum = tarenum
        tarefile = sprintf("tare%03.0f.mat",tnum);
        if isfile(strcat(taredir,'/',tarefile))
            tarefiles{tnum-1} = tarefile;

        end
        n = 1;
        
        for rnum = runnum
            % for each position
            for runpos = 1:20
                runfile = sprintf("run_%02.0f_%03.0f_%03.0f_%02.0f.mat",cnum,tnum,rnum,runpos);
                if isfile (strcat(rundir,'/',runfile))
                    runfiles{tnum-1,n} = runfile;
                    runtable{tnum-1,n} = nonDim(runfile,calfile,tarefile);
                    n = n+1;
                end
                clear runfile
            % for each power
            end
            % for each position
        end
    end
end

%%
[row,col] = size(runtable);

Vinf = zeros(row,col);
cl   = zeros(row,col);
cx   = zeros(row,col);
cm   = zeros(row,col);
rho  = zeros(row,col);
Re   = zeros(row,col);

clfig = figure('Visible','off'); hold on;
cxfig = figure('Visible','off'); hold on;
cmfig = figure('Visible','off'); hold on;
for ri = 1:row
    for ci = 1:col
        cl(ri,ci) = runtable{ri,ci}.cl_average;
        cx(ri,ci) = runtable{ri,ci}.cx_average;
        cm(ri,ci) = runtable{ri,ci}.cm_average;
        rho(ri,ci) = runtable{ri,ci}.rho;
        Vinf(ri,ci) = runtable{ri,ci}.Vinf_measured;
        
        Re(ri,ci) = (rho(ri,ci) * Vinf(ri,ci)*l)/mu;
    end
    set(0,'CurrentFigure',clfig)
    plot(cl(ri,:),Re(ri,:))
    
    set(0,'CurrentFigure',cxfig)
    plot(cx(ri,:),Re(ri,:))   
    
    set(0,'CurrentFigure',cmfig)
    plot(cm(ri,:),Re(ri,:))   
end

clfig.Visible = 'on';
cxfig.Visible = 'on';
cmfig.Visible = 'on';

% =========================================
% unblown cases
% for each run number
%     
%     load data via nondimensionalize
%     
%     Re = rho*l*V/mu;
%     
%     extract rho, Vinf, cl, cx, cm
%     
%     plot cl
%         Vinf, Re on vertical axes
%         
%     title('Re-Cl');
%     xlabel('cl');
%     
%     plot cx 
%         Vinf, Re on vertical axes
%     plot cm
%         Vinf, Re on vertical axes
%     
% end
% 
% % ================================
% % blown cases?