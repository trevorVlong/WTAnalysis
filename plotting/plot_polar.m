function plot_polar(slots,Re,dF, plot_vars)
    %Make all the plots for a specified dataset
    %dF - Flap deflection
    %slots - either 'single'/'solid' or 'slotted'
    %Re - 'H' or 'L'
    addpath('../structures/performance/')
    eval(sprintf('load %s_%s_%i.mat',slots, Re,dF));
    
    hd_c = 0.36;
    %Calculate CDadd
    [NA, NC] = size(dcj);
    cda  = zeros(NA, NC);
    tcs  = zeros(NA, NC);
    for i = 1:NA
        for j = 1:NC
            dcj_t = max(dcj(i,j),0);
            [~, ~, ~, dCJ, Tc, Vj_Vi] = get_jet_coeffs_dCJ(dcj_t, hd_c);
            tcs(i,j) = Tc;
            cda(i,j) = cx2D(i,j) + Tc;
        end
    end
    
    
    %set up colormap
%     map = jet(floor(NC));
    colormap('jet');
    cbt = '$\Delta C_J$';
    %CL-CX
    figure(1)
    hold on
    for i = 1:NA
        scatter(cx2D(i,:), cl2D(i,:), plot_vars.msize, dcj(i,:), plot_vars.marker)
    end
    hcb = colorbar();
    hcb.Title.String = cbt
    hcb.Title.Interpreter = 'latex'
    xlabel('$c_x$', 'interpreter', 'latex')
    ylabel('$c_\ell$', 'interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 14, 'FontName', 'CMU Serif')
    caxis([0 7])
    
    %CL-CM
    figure(2)
    colormap('jet');

    hold on
    for i = 1:NA
        scatter(cm2D(i,:), cl2D(i,:), plot_vars.msize, dcj(i,:), plot_vars.marker)
    end
    xlabel('$c_m$', 'interpreter', 'latex')
    ylabel('$c_\ell$', 'interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 14, 'FontName', 'CMU Serif')
    hcb = colorbar();
    hcb.Title.String = cbt
    hcb.Title.Interpreter = 'latex'
    caxis([0 7])
    %CL-CDadd
    figure(3)
    colormap('jet');
    hold on
    for i = 1:NA
        scatter(cda(i,:), cl2D(i,:), plot_vars.msize, dcj(i,:), plot_vars.marker)
    end
    xlabel('$c_{d_\mathrm{add}}$', 'interpreter', 'latex')
    ylabel('$c_\ell$', 'interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 14, 'FontName', 'CMU Serif')
    hcb = colorbar();
    hcb.Title.String = cbt
    hcb.Title.Interpreter = 'latex'
    caxis([0 7])
    %CL-a
    figure(4)
    colormap('jet');
    hold on
    for i = 1:NA
        scatter(alfa(i,:), cl2D(i,:), plot_vars.msize, dcj(i,:), plot_vars.marker)
    end
    xlabel('$\alpha$', 'interpreter', 'latex')
    ylabel('$c_\ell$', 'interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 14, 'FontName', 'CMU Serif')
    hcb = colorbar();
    hcb.Title.String = cbt
    hcb.Title.Interpreter = 'latex'
    caxis([0 7])
    %CX-a
    figure(5)
    colormap('jet');
    hold on
    for i = 1:NA
        scatter(alfa(i,:), cx2D(i,:), plot_vars.msize, dcj(i,:), plot_vars.marker)
    end
    xlabel('$\alpha$', 'interpreter', 'latex')
    ylabel('$c_x$', 'interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 14, 'FontName', 'CMU Serif')
    hcb = colorbar();
    hcb.Title.String = cbt
    hcb.Title.Interpreter = 'latex'
    caxis([0 7])
    %CM-a
    figure(6)
    colormap('jet');
    hold on
    for i = 1:NA
        scatter(alfa(i,:), cm2D(i,:), plot_vars.msize, dcj(i,:), plot_vars.marker)
    end
    xlabel('$\alpha$', 'interpreter', 'latex')
    ylabel('$c_m$', 'interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 14, 'FontName', 'CMU Serif')
    hcb = colorbar();
    hcb.Title.String = cbt;
    hcb.Title.Interpreter = 'latex';
    caxis([0 7])
    %CDadd-a
    figure(7)
    colormap('jet');
    hold on
    for i = 1:NA
        scatter(alfa(i,:), cda(i,:), plot_vars.msize, dcj(i,:), plot_vars.marker)
    end
    xlabel('$\alpha$', 'interpreter', 'latex')
    ylabel('$c_{d_\mathrm{add}}$', 'interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 14, 'FontName', 'CMU Serif')
    hcb = colorbar();
    hcb.Title.String = cbt;
    hcb.Title.Interpreter = 'latex';
    caxis([0 7])
end