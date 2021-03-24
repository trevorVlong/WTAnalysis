
function [] = setup_collection()
    %------------------------------------------------------
    % wing geometry  
    wing_geom = table();
    wing_geom.c_wing = 9 * 0.0254; % meters
    wing_geom.b_wing = 24* 0.0254; % meters
    wing_geom.prop_diam = 5 * 0.0254; % meters
    wing_geom.flap_chord = 3*0.0254; % meters
    wing_geom.R_tip = 0.0635; %meters (prop radius)
    wing_geom.r_hub = 0.014; %meters  (hub radius)

    %------------------------------------------------------
    % submodel choice
    models.motor = 'F40II-4in';
    % other options include "qprop", "F40II-4in",...
    
    %------------------------------------------------------    
    % dropbox path
    %examples
    %mac
    paths.dropboxpath = '/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/'; %adjust for you machine
    paths.runpath = strcat(paths.dropboxpath,'spring21/PerformanceRuns/runs');
    paths.tarepath = strcat(paths.dropboxpath,'spring21/PerformanceRuns/tares');
    paths.calpath = strcat(paths.dropboxpath,'spring21/PerformanceRuns/calibrations');
    paths.motorfiles = strcat(paths.dropboxpath,'/spring21/motormap/maps');
    paths.interppath = strcat(paths.dropboxpath,'/spring21/OldData');
    
    
    paths.AnalysisCode = '/Users/trevorlong/Desktop/Masters/tunnel/eSTOL';
    paths.plotting     = strcat(paths.AnalysisCode,'/plotting');
    paths.Analysis     = strcat(paths.AnalysisCode,'/Analysis');
    paths.wallcorrections = strcat(paths.Analysis,'/WallCorrections');
    paths.MotorModel      = strcat(paths.Analysis,'/MotorModel');
    paths.datamodels  = strcat(paths.Analysis,'/ReducedData');
    paths.toolboxpath = "\Users\longt\Documents\mit-git\toolbox\toolbox\aawind\aawind";
    paths.functionpath = "\Users\longt\Documents\mit-git\3by2\WTrun";
    %windows
    %dropboxpath = 'C:/users/longt/Dropbox (MIT)/Tunnel_Data';

    plotops.debug = 0; % turns debug plots on
    plotops.save = 0;  % saves
    plotops.pplot  = 0;
    plotops.close  = 1;
    plotops.savetype = 'epsc'; %typical formats are espc (Linux/Mac) or svg (Windows)

    %==========================================================================
    % plot types: 1 = on 0 = off
    plotops.contour = 1;
    plotops.VJspan_avg  = 0;
    plotops.VJplot  = 0;
    plotops.visible = "on" ; % "on" or "off"

    %------------------------------------------------------    
    % Plot labels and formatting

    plotops.font = "Avenir";
    plotops.fontstyle = "bold";
    plotops.labelFS = 12;
    plotops.titleFS = 12;
    plotops.axesFS  = 14;
    plotops.cbarFS  = 30;

    % contour plot 
    plotops.VJ_rng = 0:0.7:3.5; 
    plotops.levels = length(plotops.VJ_rng); 
    plotops.cLineWidth = 2;
    plotops.showtext = 'on';
    plotops.contourLineStyle = '-';
    plotops.colormap = 'turbo';

    % VJ/V plots
    plotops.VJlim = 6;


    % line plots
    plotops.linestyle = {'--','-.','-'};

    % scatter plots
    plotops.markerstyle = {'o','^'};


    

    %------------------------------------------------------
    % save as options.mat
    save('analysis_options.mat',...
         'models',...
         'wing_geom',...
         'paths',...
         'plotops');
end
    
    
