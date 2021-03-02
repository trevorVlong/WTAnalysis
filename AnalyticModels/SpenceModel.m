function [valtable] = SpenceModel(dcj,alfa,df,sol)
%Uses model as given by Spence (reference) to return cl(dcj,alfa,df),
%cm(dcj,alfa,df) and their associated derivatives in table format. Best way
%to use is to hand a vector of dcj,alfa,dj in column vectors of length N
%
% sol tells the program whether to use 3 point, 6 point, or 9 point fits to
% the analytical solution and takes the values of 3,6,9 accordingly


% default settings
    if nargin < 4
        sol = 3;
    end
% load settings and geometry
    load('options.mat','wing_geom');
    
    
    % convert dcj to cj
    c = wing_geom.c_wing;
    R = wing_geom.R_tip;
    rh = wing_geom.r_hub;
    b  = wing_geom.b_wing; 
    
    h = c*pi*((R/c)^2-(rh/c)^2)*4/(b/c);
    cj = dcj + 2*h/c;
    % recall dcj = 
    
    
    
    
        valtable = table();
    switch sol
        case  3
    % coeficients for 3 point solution
        A0 = 1/(4*pi) * (3.54 * cj.^(1/2) +...
                         0.325* cj + ...
                         0.156* cj.^(3/2));
                     
        B0 = 1/(4*pi) * (1.152* cj.^(1/2) +...
                         1.106* cj + ...
                         0.051* cj.^(3/2));  
                     
        % convert angles to radians
        alfa = alfa * pi/180;
        df   = df * pi/180;
        
        valtable.cl = 4*pi*A0*df + 2*pi*(1+2*B0)*alfa;
        valtable.dclda = 2*pi*(1+2*B0);
        valtable.dclddf= 4*pi*A0;
        
        case  6
            disp("I'll get around to this one");
            
        case  9
            disp("I'll get around to this one");
            
        case sol ~= 6 && sol ~= 3 && sol ~=9
            disp('choose a valid fit: 3,6,9');
    end

end

