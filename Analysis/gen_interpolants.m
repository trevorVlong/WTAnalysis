%Add cl-cm, cl-cx interpolants to existing 'single' and 'slotted' matfiles


close all
clc
addpath('/Users/trevorlong/Dropbox (MIT)/Tunnel_Data/sum_fall_20/results/structures/performance');

dF = [20, 40, 50, 55, 60];
Re = ['H', 'L'];
sl = ["single", "slotted"];
hd_c = 0.36;

for i = 1:5
    for j = 1:2 %Re
        for k = 1:2 %slots
            
            file = sprintf('%s_%s_%i.mat', sl(k), Re(j), dF(i));
            eval(sprintf('load %s', file))
            
            [NA, NC] = size(dcj);

            [cda, tc] = get_cda(dcj, cx2D, hd_c);
            
            dcj_s = reshape(dcj, 1, NA*NC);
            cl_s = reshape(cl2D, 1, NA*NC);
            cm_s = reshape(cm2D, 1, NA*NC);
            cx_s = reshape(cx2D, 1, NA*NC);
            cd_s = reshape(cda, 1, NA*NC);

            Vclcm = scatteredInterpolant(cl_s', dcj_s', cm_s');
            Vclcx = scatteredInterpolant(cl_s', dcj_s', cx_s');
            Vclcda = scatteredInterpolant(cl_s', dcj_s', cd_s');
            
            save(file, 'alfa', 'cl2D', 'cm2D', 'cx2D', 'dcj', 'flap_ang',...
                'Vcl', 'Vcm', 'Vcx',...
                'Vclcm', 'Vclcx', 'Vclcda')
        end
    end
end

xlabel('AoA')
ylabel('c_l')


