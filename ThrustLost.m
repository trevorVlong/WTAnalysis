function [ce_loss] =  ThrustLost(pred_CT, CX, rpm, rho, Vinf)
% Provides brief analysis of the losses incurred by the wing geometry based
% on thrust estimates from the wind tunnel model
%
% Inputs
% Predicted_CT := the predicted average CT for the wing written as 
%                                   CT = T/(0.5(omega * rtip)^2* A_disk)
%           CX := the net force recorded on the wing, written as 
%                                   (D-T)/(q_inf * A_wing)
% 
%
% Outputs
%  loss_coefficient = (T-X)/(.5 * rho * Vinf^3 * A_wing)

    load('analysis_options.mat','wing_geom');
    
    w = pi/30 * rpm;

    A_d = pi*(wing_geom.R_tip^2 - wing_geom.r_hub^2);
    A_wing = wing_geom.b_wing*wing_geom.c_wing;

    q_inf = 0.5 * rho * Vinf^2;
    q_tip = 0.5 * (w * wing_geom.R_tip)^2;

    X = CX * q_inf*A_wing;
    T = pred_CT * q_tip * A_d;
    
    E_loss = (T-X) * Vinf;
    
    ce_loss = E_loss/ (q_inf * Vinf * A_wing);
    
    


end
