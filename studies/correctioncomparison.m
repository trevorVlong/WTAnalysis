% compare different tunnel wall corrections


load('analysis_options.mat','wing_geom','paths','statops');
addpath(paths.runpath,paths.tarepath,paths.calpath) % adjust in "options' file

load("NominalConfig.mat");


% Mark's method
cld =[];
cxd = [];
cmd = [];
clu = [];
cxu = [];
cmu = [];
alfa = [];

clfig = figure('Visible','off');
xlabel('$\alpha$','Interpreter','latex');
ylabel('$c_\ell$','Interpreter','latex');

hold on

for jj = 1:3:11

    for ii = 1:9

        rtable = runtable{ii,jj,5};
        if istable(rtable)
            % uncorrected values
            alfa(ii) = mean(rtable.alfa) + 3.7;
            clu(ii) = cl(ii,jj,5);
            cxu(ii) = cx(ii,jj,5);
            cmu(ii) = cm(ii,jj,5);
            dcju(ii) = dcj(ii,jj,5); 

            %analytic corrections

            deltas = AnalyticWall(geom,clu(ii),cmu(ii),2);
            alfad(ii) = deltas.dalfa;
            cld(ii)   = deltas.dcl;
            %Mark's corrections

            corrd = twall(rtable);

            clM(ii) = corrd.cl_average;
        
        end
    end


    plot(alfa,clu,...
        'LineWidth',2, ...
        'Color',[0,0,0]);

    plot(alfa,clu+cld,...
        'LineWidth',2, ...
        'Color',[1,0,0]);

    plot(alfa,clM,...
        'LineWidth',2, ...
        'Color',[0,1,0]);
end
clfig.Visible = 'on';

% Infinite walls with single vortex method

