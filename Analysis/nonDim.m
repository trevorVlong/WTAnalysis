function [corrected2D,uncorrected2D] = nonDimensionalize(runfile,calfile,tarefile)
% this function corrects voltage data to non-dimensional coefficient data based on 
% the tunnel conditions and atmospheric conditions based on "weatherman" the 
% in-built weather system (in progress). Returns a vector of averaged values
%
%
% Summary:
% data a 7 x [N] vector of values from the wind tunnel that has the breakout--
%  1. Lift cell         ----> Cl
%  2. Drag cell         ----> Cx
%  3. Moment cell       ----> Cm
%  4. AoA sensor        ----> degrees
%  5. RPM1              ----> dCj1
%  6. RPM2              ----> dCj2
%  7. RPM3              ----> dCj3
%  8. RPM4              ----> dCj4
% and returns the extra channel
% these all correspond to analog inputs N-1 on the labjack
%
% Mensor ?
%
%
% weather box
% statevars.StaticPressure : static pressure in kPa
% statevars.Temperature    : Temp in celcius
% RelativeHumidity          : % relative humidity
%
%
% Inputs:
% data: a matrix as described above
% calfile: pointer to calibration file
% tarefile: pointer to tare file
% 
%
%% setup =================================================================
%options and path stuff
%save('filename.mat','-append'); to do

    load('analysis_options.mat','wing_geom','paths');
    
    addpath(paths.runpath,paths.tarepath,paths.calpath) % adjust in "options' file
    addpath('MotorModel','WallCorrections');
    
    
    
%------------------------------------
% add cal, tare files
    load(calfile,'C','alfa0V')
    load(tarefile)
    load(runfile)
    
%------------------------------------    
%add tables for use later
    uncorrected2D = table;
%corrected2D   = table;  % is now added by twall


%------------------------------------    
% pull out data as vectors (to be converted to averaged values later)
    LDMvolt         = rundata.forcedata; % Lift;Drag;Moment
    alfa            = rundata.ljdata(4,:);% angle of attack
    COUNTs          = rundata.ljdata(8:11,:); % motor RPMs
    Temp            = rundata.statevars.Temperature; % kelvin
    SPressure       = rundata.statevars.StaticPressure; 
    relHumid        = rundata.statevars.RelativeHumidity;
    tunnelq         = rundata.tunnelq;
    walldata        = rundata.walldata;
    flap_ang        = rundata.flap_ang;
    
%------------------------------------
% constants
    R = 287;
    g = 9.81;
    Mensor_offset = - 0.055;
    Mensor_multiplier = 1.869; %inches of water --> Torr
    %...
    
%----------------------------------------------------------
% % wing geometry (switch to datasheet) 11 Feb 2021
%     wing_geom = table();
%     wing_geom.c_wing = 9 * 0.0254; % meters
%     wing_geom.b_wing = 24* 0.0254; % meters
%     wing_geom.prop_diam = 5 * 0.0254; % meters
%     wing_geom.flap_chord = 3*0.0254; % meters
%     wing_geom.prop_pitch = NaN ; % idk
%     wing_geom.R_tip = 0.0635; %meters (prop radius)
%     wing_geom.r_hub = 0.014; %meters  (hub radius)

    
    
    
%% corrections

%correction factors
%-------------------------------------------------------------------------
    angle_conversion = 360/5; %volts to deg (360 deg ove 5 volt of output
    q_conversion =  2.8035*0.133322*1000; % 2.8035 Volt/Torr x .133322 kPa/Torr
%    alfa0        = 2.0739; %volts
%    alfa0        = 2.0573; %volts
    alfa0        = 2.0196;  %volts
    alfa0        = mean(alfa0V);

% state variables "statistics" and averaging
%-------------------------------------------------------------------------
% get state and freestream variables
    Temp_std     = std(Temp); % C
    Temp         = mean(Temp)+272; % C
    pambient_std = std(SPressure)*100; % pascals
    pambient     = mean(SPressure)*100; %pascals
    % no statistics b/c is dual distribution of other two
    rho          = pambient / ((Temp) * R); %kg/m^3
    
% add to data table

    
% Tunnel and Wing Conditions
%-------------------------------------------------------------------------
    % Labjack Values
    tunnelq_raw = 133.322*Mensor_multiplier*(tunnelq);
    tunnelq_std = std(tunnelq_raw);
    tunnelq = double(mean(tunnelq_raw)); % corrected down 10% while I figure out how to make q more accurate 
    % for debug
    
    % as forces
    LDM   = C*(LDMvolt-tarevec(1:3))*g;% convert voltages to forces
    %alfa0 = cal_alfa; % get 0 AoA during calibration of wing
    
% add to data table  

    
    
    % Dimensionless quantities
    %-------------------------------------------------------------------------
    % convert to coefficients
    cl     = (LDM(1,:)./ (tunnelq * wing_geom.c_wing*wing_geom.b_wing));
    cx     = (LDM(2,:)./ (tunnelq * wing_geom.c_wing*wing_geom.b_wing));
    cm     = LDM(3,:)./ (tunnelq * wing_geom.c_wing^2*wing_geom.b_wing);

    
    
    
    % motor / jet model
    udCJs = [];
    cdCJs = [];
    RPMs  = [];
    VJ_V  = [];
    uCQ   = [];
    
    try
        
        for motor = 1:4
            RPMs(motor,:) = count2RPM(COUNTs(motor,:),100);
            [udCJs(motor,:),VJ_V(motor,:),CQ(motor,:)] = RPM2DCJ(RPMs(motor,:),tunnelq,rho);
        end
        
    catch e
        %warning(" An issue occured attempting to calcualte ∆cj");
        RPMs(1:4,:) = 0;
        cdCJs(1:4,:) = 0;
        VJ_V(1:4,:) = 0;
        udCJs(1:4,:) = 0;
        CQ(1:4,:) = 0;
    end
    
%=================================================================
% table of values not corrected for wall interference

% for tunnel wall corrections 
    uncorrected2D.Fz = mean(LDM(1,:));
    uncorrected2D.Fx = mean(LDM(2,:));
    
    
% important coefficients/aero stuff
    uncorrected2D.cl_average = mean(cl);
    uncorrected2D.cx_average = mean(cx);
    uncorrected2D.cm_average = mean(cm);
    uncorrected2D.dCJ  = mean(udCJs(1:4,:),'all');
    uncorrected2D.flap = flap_ang;
    uncorrected2D.alfa = (alfa-alfa0)*angle_conversion;
    
%----------------------------------
% tunnel ambient conditions
    uncorrected2D.Vinf = sqrt(2*mean(tunnelq)/rho);
    uncorrected2D.tunnelq = tunnelq;
    uncorrected2D.Temp = Temp;
    uncorrected2D.pambient = pambient;
    uncorrected2D.rho = rho;    
    
    
%----------------------------------
%diagnostic stuff    
    
 
    % %motor performance etc
    uncorrected2D.rpm1 = RPMs(1,:);
    uncorrected2D.rpm2 = RPMs(2,:);
    uncorrected2D.rpm3 = RPMs(3,:);
    uncorrected2D.rpm4 = RPMs(4,:);
    
    uncorrected2D.dCJ1 = udCJs(1,:);
    uncorrected2D.dCJ2 = udCJs(2,:);
    uncorrected2D.dCJ3 = udCJs(3,:);
    uncorrected2D.dCJ4 = udCJs(4,:);
    uncorrected2D.uCQ  = mean(CQ,1);

    
    % deviatons
    uncorrected2D.tunnelq_std = tunnelq_std;
    uncorrected2D.alfa_std = std(alfa-alfa0)*angle_conversion;
    uncorrected2D.cl_std = std(cl);
    uncorrected2D.cx_std = std(cx);
    uncorrected2D.cm_std = std(cm);    
    uncorrected2D.Temp_std = Temp_std;
    uncorrected2D.pambient_std = pambient_std;

%==========================================================================
%table of values corrected for wall interferece

%----------------------------------
% important coefficients/aero stuff
    corrected2D = twall(uncorrected2D,walldata); %includes cl,cx Vinf_corrected
    
    % calculate new ∆cj and cm with adjusted Vinf
    for motor = 1:4
        [cdCJs(motor,:),VJ_V(motor,:),CQ(motor,:)] = RPM2DCJ(RPMs(motor,:),...
                                                               corrected2D.tunnelq,...
                                                               rho);
    end
    
    corrected2D.dCJ  = mean(cdCJs(1:4,:),'all');
    corrected2D.alfa    = mean(uncorrected2D.alfa);
    corrected2D.cm_average = uncorrected2D.cm_average * uncorrected2D.Vinf^2 / corrected2D.Vinf^2;
    
    
%----------------------------------
% tunnel ambient conditions
    corrected2D.tunnelq = .5*rho*corrected2D.Vinf^2; % update with "corrected" freestream from twall
    corrected2D.Vinf_measured = uncorrected2D.Vinf;  % old Vinf for context
    corrected2D.Temp            = uncorrected2D.Temp; %in Kelvin
    corrected2D.rho             = uncorrected2D.rho; % in kg/m^3
    corrected2D.static_pressure = uncorrected2D.pambient; % in Pa
    corrected2D.flap            = flap_ang; % in degrees
    
    
%----------------------------------
%diagnostic stuff


%motor performance, coefficients
    corrected2D.rpm1 = RPMs(1,:);
    corrected2D.rpm2 = RPMs(2,:);
    corrected2D.rpm3 = RPMs(3,:);
    corrected2D.rpm4 = RPMs(4,:);
    corrected2D.dCJ1 = cdCJs(1,:);
    corrected2D.dCJ2 = cdCJs(2,:);
    corrected2D.dCJ3 = cdCJs(3,:);
    corrected2D.dCJ4 = cdCJs(4,:);
  
    corrected2D.CQ  = mean(CQ,1);
    corrected2D.VJ_V = mean(VJ_V,'all');

% deviations
    corrected2D.cl_std = uncorrected2D.cl_std;
    corrected2D.cx_std = uncorrected2D.cx_std;
    corrected2D.cm_std = uncorrected2D.cm_std; 
    corrected2D.cm_average = uncorrected2D.cm_average;


   
end