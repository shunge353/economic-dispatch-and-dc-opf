% based on the material from DTU course 31765...
%    https://arxiv.org/pdf/1811.00943.pdf
%    https://github.com/jmontalvo94/31765_Optimization
%
clear; clc; close all;
%% load the approximate matpower case file and change to a 3 bus test system
mpc = case4gs;
define_constants;

mpc.bus(4, :) = [];
mpc.branch(4, :) = [];

% generator matrix initialize
mpc.gen(1, GEN_BUS) = 1;
mpc.gen(2, GEN_BUS) = 2;
mpc.gen(:, PMIN) = [0; 0];
% mpc.gen(:, PMAX) = [200; 200];
mpc.gen(:, PMAX) = [100; 200];    %% reduce gen 1 maximum real power output to change marginal gen

% branch matrix initialize
mpc.branch(3, T_BUS) = 3;
mpc.branch(:, [BR_R, BR_X, BR_B]) = 0;
mpc.branch(:, BR_X) = [0.1; 0.3; 0.1];
% mpc.branch(:, RATE_A) = [10000; 10000; 10000];
mpc.branch(:, RATE_A) = [10000; 40; 10000];    %% reduce branch flow limit in line 1-3 to cause congestion

% bus matrix initialize
mpc.bus(:, PD) = 0;
mpc.bus(3, PD) = 150;

% build the B matrices and phase shift injections for DC power flow
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);

c = [60; 120];    %% cost terms for generators ($/MWh)

nb = size(mpc.bus, 1);    %% number of buses
ng = size(mpc.gen, 1);    %% number of generators

baseMVA = mpc.baseMVA;
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;

%% implementation of a economic dispatch problem
Pg = sdpvar(ng, 1);

Obj = c'*Pg;

Const = [];
Const = [Const, sum(Pg) - sum(bus(:, PD)) == 0];      %% system power balance equation
Const = [Const, gen(:, PMIN) <= Pg <= gen(:, PMAX)];  %% generator real power limits

diag = sdpsettings('solver', 'gurobi', 'verbose', 0);
optimize(Const, Obj, diag);

lambda = -dual(Const(1));    %% Lagrange multiplier of system power balance equation, i.e. system marginal price ($/MWh)
disp('+++++++++++++++++++++++++++ economic dispatch ++++++++++++++++++++++++++++++');
fprintf('optimal diaptch of genertor 1 is %.2f MW, genreator 2 is %.2f MW.\n', value(Pg(1)), value(Pg(2)));
fprintf('total cost of generatos is % .2f $/h.\n', value(Obj));
fprintf('system marginal price is %.2f $/MWh.\n\n', lambda);

%% implementation of a DC OPF problem
theta = sdpvar(nb, 1);    %% theta must be expressed in radians in order to approx. sin(theta) using theta ...
Pg = sdpvar(ng, 1);       %% when voltage angle difference across branches are small enough

Obj = (c*baseMVA)'*Pg;

on = find(gen(:, GEN_STATUS) > 0);    %% which generators are on?
gbus = gen(on, GEN_BUS);              %% what buses are they at?
ngen = size(on, 1);
Cg = sparse(gbus, (1:ngen)', 1, nb, ngen);  %% connection matrix, element (i,j) is 1 if gen on(j) at bus i is ON

Const = [];
Const = [Const, Bbus*theta + Pbusinj + bus(:, PD) / baseMVA + bus(:, GS) / baseMVA - Cg*Pg == 0];  %% nodal power balance equations
Const = [Const, -branch(:, RATE_A) / baseMVA <= Bf*theta + Pfinj <= branch(:, RATE_A) / baseMVA];  %% branch flow limits
Const = [Const, 0 <= theta(1) <= 0];                                                               %% reference bus angle
Const = [Const, gen(:, PMIN) / baseMVA <= Pg <= gen(:, PMAX) / baseMVA];                           %% generator real power limits

optimize(Const, Obj, diag);

lambda = dual(Const(1)) / baseMVA;    %% Lagrange multipliers of nodal power balance equations, i.e. nodal prices ($/MWh)
cong = dual(Const(2)) / baseMVA;      %% Lagrange multipliers of branch flow limits ($/MWh)
disp('+++++++++++++++++++++++++++++++ DC OPF +++++++++++++++++++++++++++++++++++++');
fprintf('optimal dispatch of generator 1 is %.2f MW, generator 2 is %.2f MW.\n', value(Pg(1)) * baseMVA, value(Pg(2)) * baseMVA);
fprintf('total cost of generators is %.2f $/h.\n', value(Obj));
fprintf('voltage angle: bus 1 is %.2f rad, bus 2 is %.2f rad, bus 3 is %.2f rad.\n', ...
    value(theta(1)), value(theta(2)), value(theta(3)));
fprintf('nodal prices: bus 1 is %.2f $/MWh, bus 2 is %.2f $/MWh, bus 3 is %.2f $/MWh.\n', lambda(1), lambda(2), lambda(3));

if ~isempty(find(cong > 0))
    disp('congestion occurs!');
end
