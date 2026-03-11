function [out] = ConcentricWireCylinderCorona(opts,name)

if ~exist("name","var")
    name = GetTimeString(datetime);
end
p = DefaultInput(opts);

% replacing default parameters with the one specified in input structure
for folder_name = fieldnames(opts)'
    p.(folder_name{1}) = opts.(folder_name{1});
end

BCel_val = [1; 0];

%-----------------Chemistry--------------------
disp("Initializing")

[e,eps0,kB,Na,me] = FundamentalConstants();
N = 7.34e27/p.T_GAS * p.P;

if p.Kin_model == "Parent"
    nr = 13;
    qs = [-1, 0, 0, 1, 1, -1]; 
    s_names = ["e", "N2", "O2", "N2+", "O2+", "O2-"];
    ms = [me, Inf, Inf, 0.028/Na-me, 0.032/Na-me, 0.032/Na+me];
    n0 = [0.001e13, 0.7884*N, 0.2116*N, 0.8e13, 0.2e13, 0.999e13];
elseif p.Kin_model == "ParentConst"
    nr = 13;
    qs = [-1, 1, 1, -1]; 
    s_names = ["e", "N2+", "O2+", "O2-"];
    ms = [me, 0.028/Na-me, 0.032/Na-me, 0.032/Na+me];
    n0 = [0.001e13, 0.8e13, 0.2e13, 0.999e13];
elseif p.Kin_model == "Kozhevnikov"
    nr = 26;
    qs = [-1, 0, 0, 1, 1, -1, 0, 0, 1];
    s_names = ["e", "N2", "O2", "N2+", "O2+", "O2-", "O", "N", "O4+"];
    ms = [me, Inf, Inf, 0.028/Na-me, 0.032/Na-me, 0.032/Na+me, Inf, Inf, 0.064/Na-me];
    n0 = [0.001e13, 0.7884*N, 0.2116*N, 0.333e13, 0.233e13, 0.999e13, 1e10, 1e10, 0.433e13];
elseif p.Kin_model == "Townsend"
    nr = 4;
    qs = [-1, 1, -1];
    s_names = ["e", "ion+", "ion-"];
    ms = [me, 0.064/Na-me, 0.032/Na+me];
    n0 = [0.001e13, 2e10, 1e10];
end

ns = numel(s_names);

v_th = sqrt(8*kB*p.T_GAS./(pi*ms));
r = zeros(2,ns);

%---------------Mesh Generation---------------
p.a = get_a_from_np_and_delta_min(p.np,p.delta_min_relative);
r_interfaces = CreateTanhMesh(p.np,p.r0,p.R,p.a,p.b);
r_nodes = 0.5 * (r_interfaces(1:end-1) + r_interfaces(2:end));
Vol = 0.5 * (r_interfaces(2:end).^2 - r_interfaces(1:end-1).^2);
r_full = [r_interfaces(1); r_nodes; r_interfaces(end)];
delta = r_full(2:end) - r_full(1:end-1);

if p.mesh == 1
    plot(r_nodes,zeros(size(r_nodes)),"LineStyle","none","Marker",".","MarkerSize",15,"Color","r")
    xlim([r_interfaces(1),r_interfaces(end)]);
    xlabel("r (m)")
    fontsize(18,"points")

    waitfor(gcf)
end

disp('Grid info:')
fprintf('number of nodes: %d\n', p.np);
fprintf('smallest cell: %.3e (m)\n', min(diff(r_interfaces)));
fprintf('biggest cell: %.3e (m)\n\n', max(diff(r_interfaces)));

%---------------Electrostatics---------------
if p.Poisson == "Lin"
    [Kelet,ExtraRho,K_e,Add_e,E2,Diag] = Electrostatics(r_interfaces);
elseif p.Poisson == "Log"
    [~,E2,Kelet,ExtraRho,K_e,Add_e,Diag] = LogEletStat(r_interfaces);
end

N2RhoEps = Create_n_2_rho_eps_matrix(p.np,ns,qs,e,eps0);
N2RhoEps_cond = Diag * N2RhoEps;
E_Laplacian = E2(:,1);

%--------------Photoionization---------------
p_O2 = 0.2095 * p.P * 760;
Qf = 0.0395/(p.P + 0.0395);

[c1,c2,c3,A1,A2,A3] = PhotoCoefficients(p_O2);
[Kph_1,Kph_2,Kph_3,~] = CreatePhotoMatrices(p.np,r_interfaces,delta,Vol,c1,c2,c3);

%---------------Initialization---------------
if isfield(p, 'InitialCondition')
    struct = load(p.InitialCondition);
    n_start = struct.N_STEADY_STATE;
else
    n_start = ones(p.np, ns) .* n0;
end

indices_boundaries = CreateIndicesBoundaries(p.np,ns);
Gamma2puLengthCurrent = CreateSatoCurrentMatrix(p.np,ns,qs,e,r_interfaces,E_Laplacian);
Gamma2GammaClose = CreateGamma2GammaClose(p.np,ns,Vol,r_interfaces);
qs_col = reshape(ones(p.np+1,1) .* sign(qs), [], 1);
[I_D,J_D,I_Mu,J_Mu,li] = IndicesAssemble(p.np,ns);
fAssembleKD = @(D) AssembleKD(D,delta,I_D,J_D,li,p.np,ns);
fAssembleKMu = @(E,mu) AssembleKMu(E,mu,I_Mu,J_Mu,li,qs_col,p.np,ns);

Evals = load("DATA\Evals.csv");
Tevals = load("DATA\Tevals.csv");
fComputeTe = @(E) ComputeTeLoki(E, Evals, Tevals, N, p.Te);
fComputePhotoSource = @(E,s) PhotoSource(E,s,r_nodes,r_interfaces,Vol,N,Qf,Kph_1,Kph_2,Kph_3,A1,A2,A3,p.Photoionization);

if p.BC == "Absorbent"
    fApplyBC = @(n_boundary,u_boundary,~,~,~) ApplyAbsorbentBC(n_boundary,u_boundary,ns);
elseif p.BC == "GorinII"
    fApplyBC = @(n_boundary,u_boundary,Te,t,E) GorinBCII(n_boundary,u_boundary,r,p.gamma_II,qs,v_th,ns,Te);
end

if p.Kin_model == "Parent"
    fComputeCoefficients = @(E, Te) ComputeCoefficientsParent(E,Te,p.T_GAS,N,r_interfaces,r_nodes,p.MuConst);
    fComputeOmega = @(n, kr, E) ComputeOmegaParent(n, kr, E, fComputePhotoSource, p.Const_Omega);
elseif p.Kin_model == "ParentConst"
    fComputeCoefficients = @(E, Te) ComputeCoefficientsParentConst(E,Te,p.T_GAS,N,r_interfaces,r_nodes,p.MuConst);
    fComputeOmega = @(n, kr, E) ComputeOmegaParentConst(n, kr, E, N, fComputePhotoSource, p.Const_Omega);
elseif p.Kin_model == "Kozhevnikov"
    E_k1_k2 = load("DATA\E_k1_k2.csv");
    fCompute_k1_k2 = @(E) Compute_k1_k2(E, E_k1_k2, N);
    fComputeCoefficients = @(E, Te) ComputeCoefficientsKozhevnikov(E,Te,p.T_GAS,fCompute_k1_k2,N,r_interfaces,r_nodes,p.MuConst);
    fComputeOmega = @(n, kr, E) ComputeOmegaKozhevnikov(n, kr, E, fComputePhotoSource, p.Const_Omega);
elseif p.Kin_model == "Townsend"
    fComputeCoefficients = @(E, Te) ComputeCoefficientsTownsend(E,Te,p.T_GAS,N,r_interfaces,r_nodes);
    fComputeOmega = @(n, kr, E) ComputeOmegaTownsend(n, kr, E, fComputePhotoSource, p.Const_Omega);
end

if p.OutFcn == "bar"
    clear OdeProgressBar
    if max(abs(diff(p.time_instants) - mean(diff(p.time_instants)))) < 1e-10
        bar_scale = "lin";
    else
        bar_scale = "log";
    end
    output_func = @(t,y,flag) OdeProgressBar(t,y,flag,bar_scale,1);
elseif p.OutFcn == "none"
    output_func = [];
end

disp("Initialization finished")
%-------------------- Run --------------------
disp("Starting simulation")
ode_opts = odeset("OutputFcn",output_func, "RelTol",p.RelTol);
ode_opts.Mass = sparse(1:p.np*ns,1:p.np*ns,ones(1,p.np*ns),p.np*(ns+1),p.np*(ns+1));
ode_opts.MassSingular = "yes";
ode_opts.JPattern = CreateJacobianSparsityPattern(p.np,ns,qs);

phi_coherent = Kelet \ (N2RhoEps_cond * n_start(:) + ExtraRho*(BCel_val*p.v(0)));
OdeFuncDD = @(t,y) DaeFunc(t,y, ...
                           Kelet,ExtraRho,K_e,Add_e,p.np,ns,p.v,...
                           fAssembleKD,fAssembleKMu,fApplyBC,fComputeTe,fComputeCoefficients,fComputeOmega,fComputePhotoSource,...
                           Gamma2GammaClose,N2RhoEps_cond,BCel_val,indices_boundaries);
start_time = tic();
y0 = [n_start(:); phi_coherent];
ode_opts.InitialSlope = OdeFuncDD(0,y0);
[t_out,y_out,stats_out] = ode15s(OdeFuncDD,p.time_instants,y0,ode_opts);
wct = toc(start_time);
disp("Simulation finished")

%--------------Post Processing--------------
t_out = t_out';
y_out = y_out';
nt = numel(t_out);
n_out = y_out(1:p.np*ns,:);
phi_out = y_out(p.np*ns+1:end,:);
charge_density = (N2RhoEps*n_out)*eps0;
BCel_val_out = BCel_val*p.v(t_out);
E_out = K_e * phi_out + Add_e * BCel_val_out;
E_rid = (E_out/N)*1e21;
Te = ComputeTeLoki(E_out, Evals, Tevals, N, p.Te);
Sph_out = zeros(p.np,nt);
Omega_out = zeros(p.np*ns, nt);
Gamma_out = zeros((p.np+1)*ns, nt);
RR_out = zeros(p.np*nr,nt);

if p.Kin_model == "Townsend"
    KR_out = zeros(p.np*(nr-1),nt);
else
    KR_out = zeros(p.np*nr,nt);
end

for it=1:nt
    [~,Gamma,Omega,rr,Sph,kr] = OdeFuncDD(t_out(it), y_out(:,it));
    Gamma_out(:,it) = Gamma;
    Omega_out(:,it) = Omega;
    Sph_out(:,it) = Sph;
    RR_out(:,it) = rr;
    KR_out(:,it) = kr;
end

RR_out = permute(reshape(RR_out.', nt, p.np, nr), [2, 1, 3]);

if p.Kin_model == "Townsend"
    KR_out = permute(reshape(KR_out.', nt, p.np, nr-1), [2, 1, 3]);
else
    KR_out = permute(reshape(KR_out.', nt, p.np, nr), [2, 1, 3]);
end

stats.Succesful_steps = stats_out(1); 
stats.Failed_attempts = stats_out(2);
stats.Function_evaluations = stats_out(3);
stats.Partial_derivaties = stats_out(4);
stats.LU_decompositions = stats_out(5);
stats.Solutions_of_linear_systems = stats_out(6);
stats.wall_clock_time = WCTfunc(wct);

% ------------ output structure ------------
out.t_out = t_out;
out.y_out = y_out;
out.stats = stats;
out.p = p;
out.ns = ns;
out.s_names = s_names;
out.r_nodes = r_nodes;
out.r_interfaces = r_interfaces;
out.n_out = n_out;
out.rho = charge_density;
out.applied_voltage = BCel_val_out;
out.E = E_out;
out.E_rid = E_rid;
out.Te = Te;
out.I_SATO = Gamma2puLengthCurrent*Gamma_out;
out.RR = RR_out;
out.KR = KR_out;
out.Gamma = Gamma_out;
out.Omega = Omega_out;
out.Sph = Sph_out;
out.n_steadystate = reshape(n_out(:,end),p.np,ns);

Save(out, name)

end