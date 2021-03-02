% Reads wake data 
% 8 Dec, 2020
% Trevor Long
% Reads a series of wake studies specified by their calibration, tare (run series) number, and 
% ∆cj level. The dimensions of the data are also specified by Zstep and
% Nspan, which are the body z-axis step numbers and spanwise axis for the
% model. Various plots are shown and can be saved or made visible based on
% inputs from the preamble section
close all

%==========================================================================
% Setup
% use the test matrix to fill in most of this information

makestruct = 1;
flapgeom = 'slotted';
propgeom = 'large';
mountgeom = 'high';


%--------------------------------------------------------------------------
% Select Runs
load options.mat
wing_geom.prop_diam = 4*0.0254;
save('options.mat','-append','wing_geom');

tares = [27];
calnum  = 10;
dcjlevel = 1:8;

%--------------------------------------------------------------------------
% Plot active and save options

plotops.debug = 0; % turns debug plots on
plotops.save = 0;  % saves
plotops.pplot  = 1;
plotops.close  = 1;
plotops.savetype = 'epsc'; %typical formats are espc (Linux/Mac) or svg (Windows)

%--------------------------------------------------------------------------
% plot types: 1 = on 0 = off
plotops.contour = 1;
plotops.VJspan_avg  = 1;
plotops.VJplot  = 1;
plotops.visible = "on" ; % "on" or "off"

%--------------------------------------------------------------------------
% Plot labels and formatting

plotops.font = "Cmu-Serif";
plotops.fontstyle = "bold";
plotops.labelFS = 12;
plotops.titleFS = 12;
plotops.axesFS  = 14;
plotops.cbarFS  = 20;

% contour plot 
plotops.VJ_rng = 0:0.5:3; 
plotops.levels = 7; 
plotops.cLineWidth = 2;

% VJ/V plots
plotops.VJlim = 6;

% clear previous plots if asked
if plotops.close == 1
    close all
end


%--------------------------------------------------------------------------
% lookup folders



addpath(strcat(paths.dropboxpath,"/sum_fall_20/runs"));
addpath(strcat(paths.dropboxpath,"/sum_fall_20/tares"));
addpath(strcat(paths.dropboxpath,"/sum_fall_20/calibrations/"));
addpath(strcat(paths.dropboxpath,"/sum_fall_20/results/structures"));



addpath(paths.motorfiles,paths.AnalysisCode,paths.wallcorrections,paths.MotorModel);


%--------------------------------------------------------------------------
% other important folders and paths

structpath = strcat(paths.dropboxpath,"/sum_fall_20/results/structures/wakestudies");
runpath    = strcat(paths.dropboxpath,"/sum_fall_20/runs");

%--------------------------------------------------------------------------
% pitot rake geometry and steps

Zstep    = 3; % number of steps taken in z-axis during wake study
Nspan    = 19; % number of steps taken spanwise in wake study
N = Nspan*Zstep;
dzpitot = 0.125;
zpos   = 0:dzpitot:(Zstep*30*dzpitot-dzpitot);

[xpos,ypos] = meshgrid(0:Nspan,zpos);

positions = [xpos(:),ypos(:)];
%--------------------------------------------------------------------------
% build table for performance values

runtable = table();
runtable.dzpitot = dzpitot;
runtable.zpos = zpos;
runtable.Nspan = Nspan;
runtable.dspitot = 1; 


%--------------------------------------------------------------------------
% setup data structures and filenames and save structures

calfile  = sprintf("cal%02.0f.mat",calnum);

% lists for output plots
cx_c  = zeros(length(dcjlevel),3);
dcj_c = zeros(length(dcjlevel),3);
cq_c  = zeros(length(dcjlevel),3);
cl_c  = zeros(length(dcjlevel),1);

% initialize figures
figs = table();
if plotops.pplot ==1
    figs.figdcj = figure('Visible',plotops.visible);
    figs.figpolar = figure('Visible',plotops.visible);
end

tic;
for tarenum = tares
    tarefile = sprintf("tare%03.0f.mat",tarenum);
    filebase = sprintf("wake_run_%02.0f_%03.0f",calnum,tarenum);
    breakflag = 0;
    for dcjN = dcjlevel
        
        
        %build list of files to look through
        files = [];
        for runnum = 1:N
            %each location in wake
            
            
            % check if file exists
            file = sprintf("%s_%03.0f_%02.0f.mat",filebase,runnum,dcjN);
            
            if isfile(strcat(runpath,'/',file))
                files = [files;file];

            else
                breakflag = 1;
                warning('file %s does not exist, check run series and calibration number',file)
                break
            end

        end
        
        if breakflag == 1
            breakflag = 0;
            break
        end

        % save structure for wake run
        wakedata = cell(Nspan,Zstep);

        % wake pitot setup
        ambport = 31;
        Npitot  = 30;
        
        
        %==========================================================================
        % Collect and compile data

        % for averaging, est = modelled value, meas = measured value
        alfa = [];
        Dcj_est = [];
        Vinf_meas = [];
        cl_meas = [];
        cx_meas = [];
        cm_meas = [];
        VJV_est = [];
        dcj_meas = [];
        CQ_meas  = [];
        
        % loop through each file
        for nf = 1:length(files)
            %--------------------------------------------------------------------------
            % load file and retrieve data
            f = files(nf);
            load(f);
            [corrected_data,uncorrected_data] = nonDimensionalize(f,calfile,tarefile);

            %--------------------------------------------------------------------------
            % Data Reduction and error prevention
            
            % reduce wake data to VJ/V
            wd = -mean(wake_data,2);
            wakedata{nf} = wd(1:Npitot)/wd(ambport); % Cp0

            % append values for run values
            alfa = [alfa corrected_data.alfa];
            Vinf_meas = [Vinf_meas corrected_data.Vinf_measured];
            cl_meas = [cl_meas  corrected_data.cl_average];
            cx_meas = [cx_meas  corrected_data.cx_average];
            Vinf_meas = [Vinf_meas corrected_data.Vinf];
            VJV_est= [VJV_est corrected_data.VJ_V];
            CQ_meas= [CQ_meas corrected_data.CQ];

            % logic for single motor cases
            if mean(corrected_data.dCJ) < 20
                Dcj_est = [Dcj_est  corrected_data.dCJ];
            elseif mean(corrected_data.dCJ) > 20 && mean(corrected_data.dCJ1) < 20
                Dcj_est = [Dcj_est  corrected_data.dCJ1];
            elseif mean(corrected_data.dCJ1) > 20 && mean(corrected_data.dCJ2) < 20
                Dcj_est = [Dcj_est  corrected_data.dCJ2];
            elseif mean(corrected_data.dCJ2) > 20 && mean(corrected_data.dCJ3) < 20
                Dcj_est = [Dcj_est  corrected_data.dCJ3];
            elseif mean(corrected_data.dCJ3) > 20 && mean(corrected_data.dCJ4) < 20
                Dcj_est = [Dcj_est  mean(corrected_data.dCJ4)];
            end
        end
        
        
%==========================================================================
%% Data reduction        
        % add to run values table average of all runs
        runtable.alfa_avg = mean(alfa);
        runtable.Dcj_avg = mean(Dcj_est);
        runtable.cl_avg  = mean(cl_meas);
        runtable.cx_avg  = mean(cx_meas);
        runtable.VJ_V= mean(VJV_est);
        runtable.df = corrected_data.flap;
        runtable.Vinf = mean(Vinf_meas);
        runtable.rho = mean(corrected_data.rho);
        runtable.CQ  = mean(CQ_meas);
        
        
        %--------------------------------------------------------------------------
        % run readout

        fprintf('\n\n\n');
        fprintf('Run Condition and Performance\n');
        fprintf('-----------------------------------------------\n');
        fprintf('flap angle \t V_inf \t AoA \t rho\n');
        fprintf('%2.0f \t %2.1f m/s \t %2.1f deg \t %2.1f kg/m^3\n',...
            runtable.df,runtable.Vinf,runtable.alfa_avg,runtable.rho);
        fprintf('-----------------------------------------------\n');
        fprintf('(V_J/V)_est \t cl_avg \t cx_avg \t ∆Cj\n')
        fprintf('%2.1f \t %2.1f \t %2.1f \t %2.1f\n',...
            runtable.VJ_V, runtable.cl_avg, runtable.cx_avg,runtable.Dcj_avg)
        fprintf('\n')
        
        
        %--------------------------------------------------------------------------
        % Cp0 -> Vstar / Vinf 
        % Vstar = sqrt(1-Cp0+Cp) where p = pinf is assumed so it reduces to
        % Vstar/V = sqrt(1-Cp0).
        
        Cp0 = zeros(Zstep*Npitot,Nspan);
        Cp0_norm = Cp0;
        
        for Nslice = 1:Nspan
            s = [];
            
            % assemble to test segment shape
            for Nz = 1:Zstep
                s = [s;wakedata{Nslice,Nz}];
            end
            % combine in z-direction
            Cp0(:,Nslice) = s;


            % remove Cp < 0, which is evidence of massive separation
            % and static pressure loss
            
            wlen = numel(Cp0);
            for elem = 1:wlen
                if Cp0(elem) > 1
                    Cp0_norm(elem) = 1;
                else
                    Cp0_norm(elem) = Cp0(elem);
                end
            end

        end
        
        % convert to usable format of Vstar/V
        u_V = sqrt(1-Cp0_norm);

        
        %----------------------------------------------------
        % wake force and loss integrals 
        tpintegrals = wakeIntegrals(u_V,runtable,wing_geom);
        
        
        
        %----------------------------------------------------
        %readout of wake integral values, compares balance measured, full
        %integral, and 2D average integral values

        fprintf('==================================================================\n\n');
        fprintf('Comparison of Treffitz Integrals and Measurements\n');
        fprintf('------------------------------------------------\n');
        fprintf('Fx_loadcell \t Fx_integral \t Fx_2Davg\n');
        fprintf('%3.1f N \t %3.1f N \t %3.1fN\n',...
            uncorrected_data.Fx,tpintegrals.Fx_int,tpintegrals.Fxa);
        fprintf('------------------------------------------------\n');
        fprintf('\t ∆CQ \t  ∆Cj \t    Cx\n');
        fprintf('Estimated \t %+3.1f \t %+3.1f \t %+3.1f \n',...
            runtable.CQ*cosd(runtable.alfa_avg + runtable.df),runtable.Dcj_avg,runtable.cx_avg);
        fprintf('Integrated \t %+3.1f \t %+3.1f \t %+3.1f \n',...
            tpintegrals.dCQ_int,tpintegrals.dCJ_int,tpintegrals.Cx_int);
        fprintf('2D average \t %+3.1f \t %+3.1f \t %+3.1f \n',...
            tpintegrals.dCQa,tpintegrals.dCJa,tpintegrals.Cxa);
        fprintf('------------------------------------------------\n');
        [Vj_m] = spanVar(u_V,wing_geom,runtable,plotops);
        
        
        %==========================================================================
        %% plotting
        % average run values
        % need to flip() cp vector around span axis

        if plotops.pplot == 1 

            % velocity ratio ux/V
            %----------------------------------------------------------
            if plotops.VJplot == 1
                [fuV,axuV] = plotVrat(u_V,runtable,plotops);

                if plotops.save == 1
                    savename = sprintf('uV_df%2.0f_alfa%2.0f_dcj%2.0f',...
                    runtable.df,runtable.alfa_avg,runtable.Dcj_avg);
                    saveas(fuV,savename,plotops.savetype);
                end
            end
            % WORK IN PROGRESS trying to estimate flow turning here
            % assume flow leaves flap at df
            %[fduV,axduV] = plotVrat(u_V/cosd(df),runtable,plotops);



        end



       dcj_c(dcjN,:) = [runtable.Dcj_avg, tpintegrals.dCJ_int, tpintegrals.dCJa];
       cx_c(dcjN,:)  = [runtable.cx_avg, tpintegrals.Cx_int, tpintegrals.Cxa];
       cq_c(dcjN,:)  = [runtable.CQ, tpintegrals.dCQ_int, tpintegrals.dCQa];
       cl_c(dcjN,:)  = [runtable.cl_avg];
       
       
               %==========================================================================
        % save tables
        if makestruct == 1
            structname = sprintf('wake_%s_%s_%s_%02.0f_%02.0f_%02.0f',flapgeom,mountgeom,propgeom,runtable.df,runtable.alfa_avg,dcjN);
            structfile = strcat(structpath,'/',structname);
            save(structfile,'tpintegrals','u_V','Cp0_norm',...
                'Nspan','Zstep','wing_geom','runtable','files','wd','wakedata',...
                'positions')

        end
    
    end

    
%% plot comparison across runs

    if plotops.pplot == 1
        plotWakeVsBalance(figs,dcj_c,cx_c,cl_c,runtable,plotops);
    end
       %----------------------------------------------------
        % save relevant info into data structure that can be stored

end
toc