% a script that's going to demonstrate the comparison between the old motor
% data and the new data from the 1x1


%-----------------------------------------------------------
% paths

% to spring21 data
motorDropboxPath =  '/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/spring21/motormap';

% to old data
oldmotorDropboxPath =  '/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/Analysis';


%-----------------------------------------------------------
% options

ops.plot = 0;
ops.save = 0;
ops.fignumber = 0;

%-----------------------------------------------------------
%% comparison of 5 inch prop data

% comparison window
lambdarng = [0.001,0.3];
Vinf = 9;

% pull data for both

    % old data
    load(strcat(oldmotorDropboxPath,'/','ctinterp2.mat')) % only have interpolants stored from qprop
    CT5old = ctVadv(lambdarng,Vinf*ones(size(lambdarng)));
    
    
    
    % new data
% 5 inch 4 blade
    [rf5,cf5,rd5] = motorMap(20,25,1,4,ops);   
    scatter(rd5.lambda,rd5.CT,'rx')

    

% plot old interpolant, new data scatter, and new data quadratic fit
f1 = figure(1);
hold on

plot(lambdarng,CT5old,'b');
plot(cf5);



%-----------------------------------------------------------
%% comparison of 4inch 3-bladed and 5inch 4-bladed propellers




% pull data for both
    
    % 4 inch 3blade
    [rf4,cf4,rd4] = motorMap(30,34,1,3,ops);
    
    % 5 inch 4 blade
    [rf5,cf5,rd5] = motorMap(20,24,1,4,ops);
    

  
% plot fits and raw data CT-lambda

CTlambdaPlot = figure(2);
hold on
%4in
scatter(rd4.lambda,rd4.CT,'b*');
plot(cf4,'b');

%5in
scatter(rd5.lambda,rd5.CT,'rx');
plot(cf5,'r');

%qprop
plot(lambdarng,CT5old,'g');

xlabel('$\lambda$','interpreter','latex');
ylabel('$C_T$','interpreter','latex');
grid on

% plot actual thrust vs advance ratio


thrustPlot = figure(3);
hold on
% geometric stuff
d4 = convlength(5,'in','m');
d5 = convlength(5,'in','m');
w4 = pi*rd4.RPM/30;
w5 = pi*rd5.RPM/30;


rho4 = mean(rd4.rhovec);
rho5 = mean(rd5.rhovec);

qp4 = 0.5*rho4*(d4/2*w4).^2;
qp5 = 0.5*rho5*(d5/2*w5).^2;
A4  = pi*(d4/2)^2;
A5  = pi*(d5/2)^2;

T4 = rd4.CT.*qp4*A4;
T5 = rd5.CT.*qp5*A5;
T4r= rd4.Tvec;
T5r= rd5.Tvec;

%scatter(rd4.lambda,T4,'b*');
scatter(rd4.lambda,T4r,'b*');
%scatter(rd5.lambda,T5,'rx');
scatter(rd5.lambda,T5r,'rx');




