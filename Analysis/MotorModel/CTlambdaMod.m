function CT = CTlambdaMod(lambda)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/results/Analysis');
    
    load('F40-2-5in.mat','cf') % load curve fit to data
    
    a = cf.a;
    b = cf.b;
    c = cf.c;
    
    %get thrust coefficient T(lambda)
    CT = a*lambda.^2 + b*lambda + c; 
    

end