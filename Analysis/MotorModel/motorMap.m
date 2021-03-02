function [runfiles,cf,rawdata,ctInterp] = motorMap(runNumberL,runNumberH,calnumber,tarenumber,ops)
% builds an interpolant of motor map CT-lambda curves given the lowest
% runnumber and the highest runnumber for that series

    %-----------------------------------------------------------
    % defaul plotops
    if nargin < 5
        ops.plot = 1; % show plot of data
        ops.save = 0;
        ops.fignumber = 1;
        ops.propdiam = nan;
        ops.motor = 'F40-2';
    end
    
    %-----------------------------------------------------------
    % add directories to search path
    motorDropboxPath =  '/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/spring21/motormap';
    calpath  = strcat(motorDropboxPath,'/calibrations');
    tarepath = strcat(motorDropboxPath,'/tares');
    runpath  = strcat(motorDropboxPath,'/runs');
    addpath(motorDropboxPath,calpath,tarepath,runpath);
    
    %-----------------------------------------------------------
    % create file list
    runfiles = [];
    for runNumber = runNumberL:runNumberH
       runfilebase = sprintf("run_%02.0f_%03.0f_%03.0f",...
                              calnumber, tarenumber,runNumber);
       % run through power settings add to 
       for powerlevel = 1:20
           runfile = sprintf("%s_%02.0f.mat",...
                              runfilebase,powerlevel);
           % add file if it exists
           if isfile(strcat(runpath,'/',runfile))
               runfiles = [runfiles; runfile];
           end
       end

    end
    
    %-----------------------------------------------------------
    % pull data from files
    
    N = length(runfiles);
    CTvec = zeros(N,1);
    lambdavec = zeros(N,1);
    Vinfvec   = zeros(N,1);
    RPMvec    = zeros(N,1);
    Tvec      = zeros(N,1);
    rhovec    = zeros(N,1);
    qinfpropv = zeros(N,1);
    omegavec  = zeros(N,1);
    for Nrun = 1:N
        runfile = runfiles(Nrun);
        load(strcat(runpath,'/',runfile),'CT','lambda','Vinf',...
            'RPM','motorName','dprop','Xforce','rho','qinfprop','omega');
        
        % because I always miss a factor of 2 sometimes
        if runNumberL < 100
            coeff = 2;
        else
            coeff = 1;
        end
        
        
        CTvec(Nrun) = coeff*CT; % multiply by 2 for 2x and 3x numbered runs, bug in my tunnel code
        lambdavec(Nrun) = real(lambda);
        Vinfvec(Nrun) = real(Vinf);
        RPMvec(Nrun) = RPM;
        Tvec(Nrun) = -Xforce;
        rhovec(Nrun) = rho;
        qinfpropv(Nrun) = qinfprop;
        omegavec(Nrun)  = omega;
    end
    
    
    %-----------------------------------------------------------
    % curve fit/ interpolating
    
    CTvec(isnan(CTvec)) = 0;
    idx = find(CTvec == 0);
    CTvec(idx) = [];
    lambdavec(idx) = [];
    Vinfvec(idx) = [];
    RPMvec(idx) = [];
    Tvec(idx) = [];
    rhovec(idx) = [];
    qinfpropv(idx) = [];
    omegavec(idx) = [];
    
    %interpolation
    ctInterp = scatteredInterpolant(Vinfvec,RPMvec,CTvec);
    
    % 2nd order fit
    ft = fittype('a*x^2+b*x^1+c');
    %ft = fittype('a*x+b');
    [cf,gof,fitinfo] = fit(lambdavec,CTvec,ft);
    
    
    %-----------------------------------------------------------
    % plotting
    
    if ops.plot == 1
        
        % ct-lambda plot
        figure(ops.fignumber)
        hold on
        scatter(lambdavec,CTvec,'rx')
        plot(cf,'b')
        
        title(strcat('ct-lambda curve, ',motorName))
        dprop = convlength(dprop,'m','in');
        
        legend({'raw data','quadratic fit'});
    end
    

    %-----------------------------------------------------------
    %create table and save 
    
    rawdata = table();
    rawdata.CT = CTvec;
    rawdata.RPM = RPMvec;
    rawdata.lambda = lambdavec;
    rawdata.Vinf = Vinfvec;
    rawdata.Tvec = Tvec;
    rawdata.rhovec = rhovec;
    
    
    if ops.save == 1
        %savedir = '';
        savename = sprintf('%s-%1.0fin.mat',ops.motor, ops.propdiam);
        save(savename,'rawdata','cf')
    end
end

