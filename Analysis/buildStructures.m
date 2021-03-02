%function [] = buildStructures(runnum,tarenum,calnum,struct_name,Re,interpmethod)
    % Author: Trevor Long
    % Date  : 27 Oct, 2020
    % buildStructures
    % This is a script that will construct .mat files that contain raw (or
    % corrected) data as well as a linear interpolant that can be used for the
    % modelling and plotting of WT data

    % this adds relevant paths for the code to run (do not edit for individual
    % runs)
    
    

    close all
    addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/runs/')
    addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/tares/')
    addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/calibrations/')
    addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/results');
    addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/results/figures');
    addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/results/structures');
    addpath('/Users/trevorlong/Documents/MATLAB/windtunnels/3by2/WTrun/');
    savepath = "/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/spring21/OldData/";
    
    
    %--------------------------------------------
    % input section
    struct_name = "slotted";
    Re          = "H";
    method = 'natural';
    
    fplusa    = 3;
    fplusdcj  = 1;

    tarenum =  201:209;  % tare associated with this vector of flap angles
    calnum  = [10];
    dataset = [11,19]; % choose tests for that set 
    
    Ndcj    = 25; % number of dcj settings to look for

    % set up collection
    N = length(tarenum);
    uncorrected2D = [];
    dcj = [];
    cl2D = [];
    cx2D = [];
    cm2D = [];
    alfa = [];
    runfiles = [];
    checkfiles = [];
    checkerror = [];
    checknum = [];

    % create cell entries (for storage)
    cl_cell = cell(N,1);
    cx_cell = cell(N,1);
    cm_cell = cell(N,1);
    dcj_cell = cell(N,1);
    alfa_cell = cell(N,1);


    for flapN = 1:N
        %======================================================================
        % build cl,cx,cm as functions of alfa,dcj for each flap angle
    tic;

        % specify which calibration and tare files will be used
        tarefile = sprintf("tare%03.0f",tarenum(flapN));
        for runnum = dataset(1):1:dataset(2)
            row_entry = mod(runnum,10);
            %disp(" clu  cxu   cl   cx");
            for cn = calnum
                calfile  = sprintf("cal%02.0f",calnum); 


                for Nrpm = 1:Ndcj
                    tempfile = sprintf("run_%02.0f_%03.0f_%03.0f_%02.0f.mat",calnum,tarenum(flapN),runnum,Nrpm);

                    if isfile(strcat('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/runs/',tempfile))
                        runfile  = tempfile;
                        runfiles = [runfiles; runfile];

                        %tic;
                        [corrected2D] = nonDimensionalize(runfile,calfile,tarefile);
                        %ctime = toc;
                        %fprintf('corrections time: %4.2f\n',ctime);
                        dcj(row_entry,Nrpm)  = mean(corrected2D.dCJ);
                        cl2D(row_entry,Nrpm) = corrected2D.cl_average;
                        cx2D(row_entry,Nrpm) = corrected2D.cx_average;
                        cm2D(row_entry,Nrpm) = corrected2D.cm_average; 
                        alfa(row_entry,Nrpm) = mean(corrected2D.alfa);
                        flap_ang = corrected2D.flap;

                        if alfa(Nrpm) > 40
                            fprintf("%s has strange AoA data\n",runfile);
                            checkfiles = [checkfiles, runfile];
                        end


                    else
                        if Nrpm == Ndcj
                            fprintf("%s set is not a file, check entries \n \n",tempfile)
                        end
                    end
                end
            end    
        end
            % each flap setting can reset cells here

        % create interpolant structures with interpolation METHOD specified at
        % start
        
        if ~isempty(alfa)
            Vcl = scatteredInterpolant(alfa(:),dcj(:),cl2D(:),method);
            Vcx = scatteredInterpolant(alfa(:),dcj(:),cx2D(:),method);
            Vcm = scatteredInterpolant(alfa(:),dcj(:),cm2D(:),method);
            
            
            %======================================================================
            % build dcl,dcx,dcm as functions of alfa,dcj for each flap angle
            %======================================================================
            % make approximate derivatives
            dalfa = linspace(min(alfa,[],'all'),max(alfa,[],'all'),length(alfa(:,1))+fplusa);
            ddcj = linspace(0,max(dcj,[],'all'),length(dcj(:,1))+fplusdcj);
            DclDdcj = [];
            DcxDdcj = [];
            DcmDdcj = [];
            for afa = 1:length(dalfa)
                cl = Vcl(dalfa(afa)*ones(length(ddcj),1),ddcj');
                cx = Vcx(dalfa(afa)*ones(length(ddcj),1),ddcj');
                cm = Vcm(dalfa(afa)*ones(length(ddcj),1),ddcj');
                Dcl = cl(2:end)-cl(1:end-1);
                Dcx = cx(2:end)-cx(1:end-1);
                Dcm = cm(2:end)-cm(1:end-1);
                Ddcj = ddcj(2:end)-ddcj(1:end-1);

                DclDdcj(afa,:) = Dcl./Ddcj';
                DcxDdcj(afa,:) = Dcx./Ddcj';  
                DcmDdcj(afa,:) = Dcm./Ddcj';  
            end

            filename = sprintf("%s_%s_%2.0f",struct_name,Re,flap_ang);
            filepath = strcat(savepath,filename);
            save(filepath,'Vcl','Vcx','Vcm','cl2D',...
                          'cx2D','cm2D','alfa','dcj',...
                          'flap_ang','runfiles','DclDdcj',...
                          'DcxDdcj','DcmDdcj');

        end
        
        ftime = toc;
        fprintf('flap deflection %2.0f runtime: %4.2f\n',flap_ang,ftime);
    end
%end