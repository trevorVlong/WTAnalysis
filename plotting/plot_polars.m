%%%%%%%%%%%%
%Plot polars
%%%%%%%%%%%%
clear all;
close all;
dF = [20, 40, 50, 55, 60];

%Plots the following figures for different combinations of datasets
%1: cl-cx
%2: cl-cm
%3: cl - cda
%4: cl - a
%5: cx - a
%6: cm - a
%7: cda - a
%cda = cx + ct;

%%%%%%%%%%%%%%%%%%%%%%%%
%Compare Re Data
%%%%%%%%%%%%%%%%%%%%%%%%
plot_vars.marker = 'o';
plot_vars.msize = 25;
plot_polar('solid', 'H', 40, plot_vars)
plot_vars.marker = 'x';
plot_polar('solid', 'L', 40, plot_vars)

figure(7)
p(1) = scatter(0,0,15, 'k', 'o');
p(2) = scatter(0,0,15, 'k', 'x');
legend(p, 'Re: 0.15e6', 'Re: 0.1e6')
%saveas(gcf,'RE_cmp/Cda_a_40', 'epsc');
figure(6)
p(1) = scatter(0,0,15, 'k', 'o');
p(2) = scatter(0,0,15, 'k', 'x');
legend(p, 'Re: 0.15e6', 'Re: 0.1e6')
%saveas(gcf,'RE_cmp/Cm_A_40', 'epsc');

figure(5)
p(1) = scatter(0,0,15, 'k', 'o');
p(2) = scatter(0,0,15, 'k', 'x');
legend(p, 'Re: 0.15e6', 'Re: 0.1e6')
%saveas(gcf,'RE_cmp/CX_A_40', 'epsc');

figure(4)
p(1) = scatter(0,0,15, 'k', 'o');
p(2) = scatter(0,0,15, 'k', 'x');
legend(p, 'Re: 0.15e6', 'Re: 0.1e6')
%saveas(gcf,'RE_cmp/CL_a_40', 'epsc');

figure(3)
p(1) = scatter(0,0,15, 'k', 'o');
p(2) = scatter(0,0,15, 'k', 'x');
legend(p, 'Re: 0.15e6', 'Re: 0.1e6')
%saveas(gcf,'RE_cmp/CL_CDA_40', 'epsc');

figure(2)
p(1) = scatter(0,0,15, 'k', 'o');
p(2) = scatter(0,0,15, 'k', 'x');
legend(p, 'Re: 0.15e6', 'Re: 0.1e6')
%saveas(gcf,'RE_cmp/CL_CM_40', 'epsc');

figure(1)
p(1) = scatter(0,0,15, 'k', 'o');
p(2) = scatter(0,0,15, 'k', 'x');
legend(p, 'Re: 0.15e6', 'Re: 0.1e6')
%saveas(gcf,'RE_cmp/CL_CX_40', 'epsc');

