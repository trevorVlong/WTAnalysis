function [outdata1,outdata2] = plotInterpolantPolars(handle,alfadata,dcjdata,interpolant1,interpolant2,n)

load('analysis_options.mat','plotops');
Ndcj = length(dcjdata);
Nalfa = length(alfadata);


% colormap setup
k = 1000;
cmap = turbo(k);
dcjcolor_index = floor((dcjdata/(7))*1000)+1;


    for dcjN = 1:Ndcj
        dcj = dcjdata(dcjN)*ones(size(alfadata));
        valdata1 = interpolant1(alfadata,dcj);
        valdata2 = interpolant2(alfadata,dcj);
        outdata1(dcjN,:) = valdata1;
        outdata2(dcjN,:) = valdata2;
        plot(handle,valdata1,valdata2,...
                    'Color',cmap(dcjcolor_index(dcjN),:),...
                    'LineStyle',plotops.linestyle{n},...
                    'LineWidth',plotops.linewidth,...
                    'HandleVisibility','off' )

        hold on;
    end
end