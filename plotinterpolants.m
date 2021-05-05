function [outdata] = plotinterpolants(handle,alfadata,dcjdata,interpolant,n)

load('analysis_options.mat','plotops');
Ndcj = length(dcjdata);
Nalfa = length(alfadata);


% colormap setup
k = 1000;
cmap = turbo(k);
dcjcolor_index = floor((dcjdata/(7))*1000)+1;


    for dcjN = 1:Ndcj
        dcj = dcjdata(dcjN)*ones(size(alfadata));
        valdata = interpolant(alfadata,dcj);
        outdata(dcjN,:) = valdata;
        plot(handle,alfadata,valdata,...
                    'Color',cmap(dcjcolor_index(dcjN),:),...
                    'LineStyle',plotops.linestyle{n},...
                    'LineWidth',plotops.linewidth,...
                    'HandleVisibility','off' )

        hold on;
    end
end