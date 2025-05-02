%--------------------------------------------------------------------------
% 09/30/24, K.K., initial version
% 04/16/25, K.K., update version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fits Tribolium beetle population data using the LA model.
% This model captures the dynamics of Larvae (L) and Adults (A).
% Parameters are estimated using optimization with MultiStart.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

%% Load Data
rng('default')
warning('off','all')
Data = readtable("TriboliumData_2.xlsx"); % Aggregated data 
warning('on','all')

TimeStep = 1:max(Data.Time);
Week = (TimeStep - 1) * 2;  % Convert time steps to weeks (starting at week 0)
N1 = max(Data.Time);        % Number of time steps

%% Calculate Means
Means = table;
Means.Larvae = zeros(N1,1);
Means.Pupae  = zeros(N1,1);
Means.Adult  = zeros(N1,1);

for i = TimeStep
    Means.Larvae(i) = mean(Data.Larvae(Data.Time == i));
    Means.Pupae(i)  = mean(Data.Pupae(Data.Time == i)); 
    Means.Adult(i)  = mean(Data.Adult(Data.Time == i)) + ...
                      mean(Data.Callows(Data.Time == i)); % Combine adults and callows
end

%% Calculate Standard Deviations
sd = table;
sd.Larvae = zeros(N1,1);
sd.Pupae  = zeros(N1,1);
sd.Adult  = zeros(N1,1);

for i = TimeStep
    sd.Larvae(i) = std(Data.Larvae(Data.Time == i));
    sd.Pupae(i)  = std(Data.Pupae(Data.Time == i)); 
    sd.Adult(i)  = std(Data.Adult(Data.Time == i)) + ...
                   std(Data.Callows(Data.Time == i));
end

%% Compute Confidence Intervals (approx. 95%)
CI = table;
CI.Larvae = (2.201 * sd.Larvae) / sqrt(12);
CI.Pupae  = (2.201 * sd.Pupae) / sqrt(12);
CI.Adult  = (2.201 * sd.Adult) / sqrt(12);

% Combine Larvae and Pupae for model fitting
Ldata = Means.Larvae + Means.Pupae;
Adata = Means.Adult;

%% Parameter Optimization
%    [b,   c1,    c2,  mu_a, mu_l]
lb = [1    0      0    0     0];     % Lower bounds
ub = [60   0.2    0.2  1     1];     % Upper bounds
p0 = [10   0.01   0.01 0.2   0.2];   % Initial guesses

problem = createOptimProblem('fmincon', ...
    'x0', p0, ...
    'objective', @(params)obj_fun(params, Ldata, Adata, N1), ...
    'lb', lb, 'ub', ub);

ms = MultiStart("Display", "iter");
[best_params_ms, fval] = run(ms, problem, 1000);

% Simulate model using best-fit parameters
[Lt, At, Q_LA] = LA(best_params_ms, Ldata, Adata, N1);
fprintf('\nBest-fit parameters:\n');
fprintf('  b     = %.4f\n', best_params_ms(1));
fprintf('  c1    = %.4f\n', best_params_ms(2));
fprintf('  c2    = %.4f\n', best_params_ms(3));
fprintf('  mu_a  = %.4f\n', best_params_ms(4));
fprintf('  mu_l  = %.4f\n', best_params_ms(5));

fprintf('\nEstimation error (relative L1 norm):\n');
fprintf('  Q_L   = %.4f\n', Q_LA(1));
fprintf('  Q_A   = %.4f\n\n', Q_LA(2));


%% Plotting: Model vs Data
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','none');
set(0, 'defaulttextInterpreter','none', 'defaultAxesFontSize', 18)
format long

Week = (TimeStep - 1) * 2;

figure()
t = tiledlayout(2,1);
nexttile
hold on
errorbar(Week(1:length(Ldata)), Ldata, CI.Larvae, 'bo', ...
    'MarkerSize', 8, 'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', [0.7 0.7 1], 'LineWidth', 2);
plot(Week, Lt, 'b.-', 'LineWidth', 2)
xlim([Week(1) 34])
box on
ylabel("Larvae", 'FontName', 'Times New Roman')
set(gca, 'FontName', 'Times New Roman')
title('LA')

nexttile
hold on
errorbar(Week(1:length(Adata)), Adata, CI.Adult, 'ko', ...
    'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'Color', 'k', 'LineWidth', 2);
plot(Week, At, 'k.-', 'LineWidth', 2)
xlim([Week(1) 34])
box on
ylabel("Adult", 'FontName', 'Times New Roman')
xlabel("Week", 'FontName', 'Times New Roman')
lg = legend("Data", "Model output", 'numcolumns', 2, 'FontName', 'Times New Roman');
lg.Layout.Tile = 'South';
set(gca, 'FontName', 'Times New Roman')

%% Objective Function: Error between data and model
function err = obj_fun(params, Ldata, Adata, Nmax)
    [~, ~, Q] = LA(params, Ldata, Adata, Nmax);
    err = 0.5 * Q(1) + 0.5 * Q(2);  % Weighted average error
end

%% LA Model Simulation Function
function [L, A, Q] = LA(params, Ldata, Adata, Nmax)
    % Extract parameters
    b = params(1);
    c1 = params(2); 
    c2 = params(3);
    mu_a = params(4);
    mu_l = params(5);

    L = zeros(Nmax,1);
    A = zeros(Nmax,1);

    L(1) = Ldata(1);
    A(1) = Adata(1);

    for n = 1:Nmax-1
        L(n+1) = A(n)*b*exp(-c1*A(n)) + (1 - mu_l)*L(n);
        % Alternative with L effect: 
        % L(n+1) = A(n)*b*exp(-c1*A(n)-c2*L(n)) + (1 - mu_l)*L(n);
        A(n+1) = mu_l*L(n)*exp(-c2*A(n)) + A(n)*(1 - mu_a);
    end

    % Relative error (normalized L1)
    Q = zeros(2,1);
    if Nmax == length(Ldata)
        Q(1) = sum(abs(Ldata(2:end) - L(2:end))) / sum(Ldata(2:end) + L(2:end));
        Q(2) = sum(abs(Adata(2:end) - A(2:end))) / sum(Adata(2:end) + A(2:end));
    end
end
