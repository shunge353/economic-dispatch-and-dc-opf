% based on the material from DTU course 31765...
%    https://arxiv.org/pdf/1811.00943.pdf
%    https://github.com/jmontalvo94/31765_Optimization
%
clear; clc; close all
%% load the approximate matpower case file and define some constants
mpc = swiss_dcopf_LP;
define_constants;

nb = size(mpc.bus, 1);    %% number of buses
ng = size(mpc.gen, 1);    %% number of generators

baseMVA = mpc.baseMVA;
bus = mpc.bus;
branch = mpc.branch;
gen = mpc.gen;
gencost = mpc.gencost;

% bus(:, PD) = bus(:, PD) * .7;  %% reduce system load to 70% of the nominal loading

[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);  %% build the B matrices and phase shift injections for DC power flow

%% implementation of a DC OPF problem
Pg = sdpvar(ng, 1);
theta = sdpvar(nb, 1);

Obj = (gencost(:, COST)*baseMVA)'*Pg;

on = find(gen(:, GEN_STATUS > 0));    %% which generators are on?
gbus = gen(on, GEN_BUS);              %% what buses are they at?
ngen = size(on, 1);
Cg = sparse(gbus, (1:ngen)', 1, nb, ng);  %% connection matrix, element (i,j) is 1 if gen on(j) at bus i is ON

Const = [];
Const = [Const, Bbus*theta + Pbusinj + bus(:, PD) / baseMVA + bus(:, GS) / baseMVA - Cg*Pg == 0];  %% nodal power balance equations
Const = [Const, -branch(:, RATE_A) / baseMVA <= Bf*theta + Pfinj <= branch(:, RATE_A) / baseMVA];  %% branch flow limits
Const = [Const, 0 <= theta(1) <= 0];                                                               %% reference bus angle
Const = [Const, gen(:, PMIN) / baseMVA <= Pg <= gen(:, PMAX) / baseMVA];                           %% generator real power limits

diag = sdpsettings('solver', 'gurobi', 'verbose', 0);

optimize(Const, Obj, diag);

lambda = dual(Const(1)) / baseMVA;    %% Lagrange multipliers of nodal power balance equations, i.e. nodal prices ($/MWh)
cong = dual(Const(2)) / baseMVA;      %% Lagrange multipliers of branch flow limits ($/MWh)
Pg = value(Pg) * baseMVA;             %% generators real power output (MW)
disp('++++++++++++++++++++++++++++++ DC OPF ++++++++++++++++++++++++++++++++++');
fprintf('total cost of generators is %.2f$/h.\n', value(Obj));

if ~isempty(find(cong > 0))
    disp('congestion occurs!');
end

%% DC OPF parameterized by system load and gen 5 maximum real power output
Pd = sdpvar(nb, 1);
Pg5_max = sdpvar(1,1);
Pg = sdpvar(ng, 1);
theta = sdpvar(nb, 1);

Obj = (gencost(:, COST)*baseMVA)'*Pg;

Const = [];
Const = [Const, Bbus*theta + Pbusinj + Pd / baseMVA + bus(:, GS) / baseMVA - Cg*Pg == 0];  %% nodal power balance equations
Const = [Const, -branch(:, RATE_A) / baseMVA <= Bf*theta + Pfinj <= branch(:, RATE_A) / baseMVA];  %% branch flow limits
Const = [Const, 0 <= theta(1) <= 0];                                                               %% reference bus angle
Const = [Const, gen(1:ng-1, PMIN) / baseMVA <= Pg(1:ng-1, 1) <= gen(1:ng-1, PMAX) / baseMVA];                           %% generator real power limits
Const = [Const, gen(ng, PMIN) / baseMVA <= Pg(ng, 1) <= Pg5_max / baseMVA];

diag = sdpsettings('solver', 'gurobi', 'verbose', 0);

P = optimizer(Const, Obj, diag, {Pd, Pg5_max}, Obj);

for i = .2:.1:1
    for j = 0:.1:1
        Pd = bus(:, PD)*i;
        Pg5_max = gen(5, PMAX)*j;
        sol= P({Pd, Pg5_max});
    end
end
