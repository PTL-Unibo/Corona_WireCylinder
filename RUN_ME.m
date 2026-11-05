addpath(genpath("src"))
schemes = ["3-species", "6-species", "9-species", "PTL-AirCorona"];
for s = schemes
    if ~isfolder(fullfile("test_case", s))
        mkdir(fullfile("test_case", s))
    end
end

start_time = tic();

% Default input -----------------------------------------------------------
opts.r0 = 300e-6;
opts.R = 3.75e-2;
opts.np = 400;
opts.delta_min_relative = 5e-5;
opts.time_instants = logspace(-12,-2,200);
opts.BC = 'GorinII';
opts.Photoionization = 1;
opts.P = 1;

% Applied voltages for each kinetic scheme --------------------------------
% (same order as schemes)
V = {(15:0.5:25)*1e3, ...  % 3-species
     (10:0.5:25)*1e3, ...  % 6-species
     (10:0.5:25)*1e3, ...  % 9-species
     (10:0.5:25)*1e3};     % PTL-AirCorona

% I-V curve using the 3-species kinetic scheme ----------------------------
opts.Kin_model = 'Townsend';
v = V{1};
for i = 1:length(v)
    opts.v = @(t) LinRamp(t, 1e-4, 10e3, v(i));
    out = ConcentricWireCylinderCorona(opts,"test_case/3-species/"+v(i));
end

% I-V curve using the 6-species kinetic scheme ----------------------------
opts.Kin_model = 'Parent';
v = V{2};
for i = 1:length(v)
    opts.v = @(t) LinRamp(t, 1e-4, v(1), v(i));
    out = ConcentricWireCylinderCorona(opts,"test_case/6-species/"+v(i));
end

% I-V curve using the 9-species kinetic scheme ----------------------------
opts.Kin_model = 'Kozhevnikov';
v = V{3};
for i = 1:length(v)
    opts.v = @(t) LinRamp(t, 1e-4, v(1), v(i));
    out = ConcentricWireCylinderCorona(opts,"test_case/9-species/"+v(i));
end

% I-V curve using the PTL-AirCorona kinetic scheme ------------------------
% (use opts.Kin_model = 'PTL-AirCorona_reduced' for the reduced scheme)
opts.Kin_model = 'PTL-AirCorona';
v = V{4};
for i = 1:length(v)
    opts.v = @(t) LinRamp(t, 1e-4, v(1), v(i));
    out = ConcentricWireCylinderCorona(opts,"test_case/PTL-AirCorona/"+v(i));
end

wct = toc(start_time);
fprintf("%s\n", "Wall clock time: ", WCTfunc(wct))

% Plot simulation results -------------------------------------------------
% each scheme is plotted over its own voltage range
exp_data = load(fullfile("DATA","300um_1_atm.csv"));
colors = [230, 25, 25; 25, 90, 230; 80, 80, 80; 255, 140, 0]/255;

figure; hold on
for k = 1:numel(schemes)
    v = V{k};
    I = zeros(size(v));
    for i = 1:length(v)
        out = load("test_case/"+schemes(k)+"/"+v(i)+".mat", "I_SATO");
        I(i) = out.I_SATO(end);
    end
    plot(v/1e3, I*1e3, "LineWidth", 2, "Color", colors(k,:));
end
plot(exp_data(:,2)/1e3, (exp_data(:,1)/0.2)*1e3, "Marker", ".", "MarkerSize", 20, ...
     "LineStyle", "none", "Color", [30, 30, 30]/255);

ax = gca;
ax.YLim = [0,3.5];
ax.YTick = 0:0.5:3.5;
ax.XLim = [10,26];
ax.XTick = 10:26;
xlabel("V (kV)")
ylabel("I (mA/m)")

legend("3-species model","6-species model","9-species model","PTL-AirCorona","Exp.","Location","northwest")
fontsize(20, "points")
