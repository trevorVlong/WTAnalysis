%plot grids of all polars
flaps = [20,40,50,55,60];
types = ["solid_mod"];
plotvars.names = ["Nominal Mount"];
Re    ="L";
plotvars.save = ""; % "save" to save plots
plotvars.linestyles = ['-','-.'];

%% polar
plotvars.fh = figure(1);
plotvars.Msize = 25;
plotvars.xlabel = "$c_x$";
plotvars.xlim = [-5 1];
plotvars.ylabel = "$c_l$";
plotvars.ylim = [-1 11];
plotvars.suptitle = "X-force Polar, V_{\infty} = 23mph";
plotvars.fsize  = 14;
plotvars.type = "polar";
plotvars.plotname = "polar_comparison";
plotData(flaps,types,Re,plotvars)

%% cl-alpha
plotvars.fh = figure(2);
plotvars.Msize = 25;
plotvars.xlabel = "$\alpha$";
plotvars.xlim = [-5 25];
plotvars.ylabel = "$c_l$";
plotvars.ylim = [-1 11];
plotvars.suptitle = "Lift Coefficient, V_{\infty} = 23mph";
plotvars.type = "cla";
plotvars.fsize  = 14;
plotvars.plotname = "cla_comparison";
plotData(flaps,types,Re,plotvars)

%% cx - alpha
plotvars.fh = figure(3);
plotvars.Msize = 25;
plotvars.xlabel = "$\alpha$";
plotvars.xlim = [-5 25];
plotvars.ylabel = "$c_x$";
plotvars.ylim = [-6 1];
plotvars.suptitle = "X-force Coefficient, V_{\infty} = 23mph";
plotvars.type = "cxa";
plotvars.fsize  = 14;
plotvars.plotname = "cxa_comparison";
plotData(flaps,types,Re,plotvars)

%% cm - alpha
Re = "L";
plotvars.fh = figure(4);
plotvars.Msize = 25;
plotvars.xlabel = "$\alpha$";
plotvars.xlim = [-5 25];
plotvars.ylabel = "$c_m$";
plotvars.ylim = [-1.2 1];
plotvars.suptitle = "Moment Coefficient, V_{\infty} = 15mph";
plotvars.type = "cma";
plotvars.fsize  = 14;
plotvars.plotname = "cma_comparison";
plotData(flaps,types,Re,plotvars)


