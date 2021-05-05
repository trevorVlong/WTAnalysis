% load file
%filename = '/Users/grace/Downloads/';
%load(filename);






%%
colormap jet
sz = size(perftable.cx);
cx = perftable.cx;
cx(cx == 0) = nan;
alfa = perftable.alfa;
alfa(alfa==0) = nan;
set(0,'DefaultAxesColorOrder',jet(11))
dcj = perftable.dcj;
dcj(dcj > 10) = nan;

for row = 2:sz(1)
    
    cx_row = perftable.cx(row,1,1);
    cx_row = cx_row(:);
    min_cx = min(cx_row);
    delta_cx = -(min_cx - perftable.cx(row,:,:)) ./ abs(min_cx);
    % get rid of big boi vals
    delta_cx(abs(delta_cx) > 30) = 0;
    alfa(alfa>30) = nan;
    
    for graph_num = 1:9
        df = perftable.df(1,1,graph_num);
        subplot(3,3,graph_num)
        %plot(perftable.alfa(:,:, graph_num), delta_cx(:,:, graph_num))
        perftable.alfa(perftable.alfa > 30) = 0;
        scatter(alfa(row,:, graph_num), delta_cx(1,:,graph_num), 20,dcj(row,:,graph_num),'filled')%, 20, perftable.dcj(1,:,graph_num), 'fill')
        xlim([-15,30])
        ylim([-.1 2])
        hold on
        title(sprintf('$\\delta_f$ = %2.0f',df),'Interpreter','latex');
        xlabel('Angle of Attack');
        ylabel('loss fraction');
        colorbar
    end
end