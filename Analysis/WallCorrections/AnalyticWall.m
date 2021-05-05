function [deltas] = AnalyticWall(geom,cl,cm,order)
%Applies results from Analytic corrections for a ~small~ 2D wing in a 
%wind tunnel that are taken from
%https://apps.dtic.mil/sti/pdfs/ADA356695.pdf 
%
% Assumptions
%
% - 2D flow
% - c_airfoil << H (smaller is relative, there is no clear ratio of where
% this breaks down
% - infinite walls (no unloading of airfoil a. la. open jet)
%
% inputs
% cl: recorded cl Nx1 vector
% cm: recorded cm Nx1 vector
% geom: struct containing wind tunnel and model dimensions uses geom.H
% (scalar)
% order: [2,4] correction order for (c/h), 2nd and 4th order solutions can be
% used

if nargin <4
    order = 2;
end

% setup
deltas = struct;
h      = geom.H; % in m
c      = geom.chord; % in m

% mach number correction (1-M^2)^0.5, set to M = 0 for assumed
% incompressible flow
beta   = 1; 


%% corrections

switch order
    
    % second order
    % innacurate for c > 0.4*beta*H
    case 2
        
        %alfa correction ∆a (radians)
        deltas.dalfa = (pi*c^2)/(96*beta*h^2)*(cl + 4*cm);
        
        % cl correction ∆c_l
        deltas.dcl = - pi^2/48*(c/(beta*h))^2 * cl;
        
        % cm correction
        deltas.dcm = pi^2/192 * (c/(beta*h))^2 * cl;
        
        % second order
    case 4
        
        %alfa correction ∆a (radians)
        deltas.dalfa = (pi*c^2)/(96*beta*h^2)*(cl + 4*cm) ...
            - (7 * pi^3*c^4 * cl)/ (30720*beta^3*h^4);
        % convert to degrees
        deltas.dalfa = 180/pi * deltas.dalfa;
        
        % cl correction ∆c_l
        deltas.dcl = cl*(-( pi^2/48 * (c/(beta*h))^2)...
                     +  (7*pi^4/3072 * (c/(beta*h)^4))  );
        
        % cm correction
        deltas.dcm = cl*( ( pi^2/192 * (c/(beta*h))^2)...
                     -  (7*pi^4/15360 * (c/(beta*h)^4))  );
    otherwise 
        warning('Not a valid correction order, only orders [2,4] are supported');
        
end

