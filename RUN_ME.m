addpath(genpath("src"))
mkdir('Test_case')
mkdir('Test_case','3-species')
mkdir('Test_case','6-species')
mkdir('Test_case','9-species')

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

% I-V curve using the 3-species kinetic scheme ----------------------------
opts.Kin_model = 'Townsend';
v = (10:0.5:25)*1e3;
for i = 1:length(v)
    opts.v = @(t) LinRamp(t, 1e-4, v(1), v(i));
    out = ConcentricWireCylinderCorona(opts,"Test_case/3-species/"+v(i));
end

% I-V curve using the 6-species kinetic scheme ----------------------------
opts.Kin_model = 'Parent';
v = (10:0.5:25)*1e3;
for i = 1:length(v)
    opts.v = @(t) LinRamp(t, 1e-4, v(1), v(i));
    out = ConcentricWireCylinderCorona(opts,"Test_case/6-species/"+v(i));
end

% I-V curve using the 9-species kinetic scheme ----------------------------
opts.Kin_model = 'Kozhevnikov';

v = (10:0.5:25)*1e3;
for i = 1:length(v)
    opts.v = @(t) LinRamp(t, 1e-4, v(1), v(i));
    out = ConcentricWireCylinderCorona(opts,"Test_case/9-species/"+v(i));
end

wct = toc(start_time);
fprintf("%s\n", "Wall clock time: ", WCTfunc(wct))

% Plot simulation results -------------------------------------------------
exp = load("DATA/300um_1_atm.csv");
I3 = [];
I6 = [];
I9 = [];
v = (10:0.5:25)*1e3;

for i = 1:length(v)
    out3 = load("Test_case/3-species/"+v(i));
    out6 = load("Test_case/6-species/"+v(i));
    out9 = load("Test_case/9-species/"+v(i));
    I3 = [I3, out3.I_SATO(end)];
    I6 = [I6, out6.I_SATO(end)];
    I9 = [I9, out9.I_SATO(end)];
end

figure; hold on
p = plot(v/1e3,I3*1e3,...
         v/1e3,I6*1e3,...
         v/1e3,I9*1e3,"LineWidth",2);
exp = plot(exp(:,2)/1e3,(exp(:,1)/0.2)*1e3,"Marker",".","MarkerSize",20,"LineStyle","none");
p(1).Color = [230, 25, 25]/255;
p(2).Color = [ 25, 90, 230]/255;
p(3).Color = [255, 140, 0]/255;
exp.Color = [30, 30, 30]/255;

ax = gca;
ax.YLim = [0,3.5];
ax.YTick = 0:0.5:3.5;
ax.XLim = [10,26];
ax.XTick = 10:26;
xlabel("V (kV)")
ylabel("I (mA/m)")

legend("3-species model","6-species model","9-species model","Exp.","Location","northwest")
fontsize(20, "points")