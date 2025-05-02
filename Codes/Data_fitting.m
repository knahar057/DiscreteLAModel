%--------------------------------------------------------------------------
% 09/30/24, K.K., initial version
% 04/16/25, K.K., update version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fits Alme data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

%% Load Data

rng('default')                    % For reproducibility
warning('off','all')              % Suppress warnings during data loading
Data = readtable("TriboliumData_2.xlsx");
warning('on','all')               % Re-enable warnings

TimeStep = 1:max(Data.Time);      % Define time steps based on data
Week = (TimeStep - 1) * 2;        % Convert to weeks, starting at week 0
N1 = max(Data.Time);              % Total number of time steps

%% Compute Means

Means = table;                    % Preallocate table for means
Means.Larvae = zeros(N1,1);
Means.Pupae = zeros(N1,1);
Means.Adult = zeros(N1,1);

for i = TimeStep
    Means.Larvae(i) = mean(Data.Larvae(Data.Time == i));
    Means.Pupae(i) = mean(Data.Pupae(Data.Time == i), "omitnan");  % Omit missing value at step 5
    Means.Adult(i)  = mean(Data.Adult(Data.Time == i)) + ...
                      mean(Data.Callows(Data.Time == i));
end

%% Compute Standard Deviations

sd = table;                       % Preallocate table for standard deviations
sd.Larvae = zeros(N1,1);
sd.Pupae = zeros(N1,1);
sd.Adult = zeros(N1,1);

for i = TimeStep
    sd.Larvae(i) = std(Data.Larvae(Data.Time == i));
    sd.Pupae(i) = std(Data.Pupae(Data.Time == i), "omitnan");      % Omit missing value at step 5
    sd.Adult(i)  = std(Data.Adult(Data.Time == i)) + ...
                   std(Data.Callows(Data.Time == i));
end

%% Compute Confidence Intervals (95%)

CI = table;                       % Preallocate table for confidence intervals
CI.Larvae = (2.201 * sd.Larvae) / sqrt(12);
CI.Pupae  = (2.201 * sd.Pupae)  / sqrt(12);
CI.Adult  = (2.201 * sd.Adult)  / sqrt(12);

% Store data vectors for plotting
Ldata = Means.Larvae;
Pdata = Means.Pupae;
Adata = Means.Adult;

%% Plot Data

% Set default interpreters and font size
set(groot, 'defaultAxesTickLabelInterpreter', 'none');
set(groot, 'defaultLegendInterpreter', 'none');
set(0, 'defaulttextInterpreter', 'none', 'defaultAxesFontSize', 18);
format long

% Regenerate week vector in case needed again
Week = (TimeStep - 1) * 2;

figure()
title("Data for PN10", 'Interpreter', 'latex', 'FontSize', 16);
hold on

% Plot Larvae
errorbar(Week(1:length(Ldata)), Ldata, CI.Larvae, 'bo-', ...
    'MarkerSize', 8, 'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', [0.7 0.7 1], 'LineWidth', 2);

% Plot Pupae
errorbar(Week(1:length(Ldata)), Pdata, CI.Pupae, 'ro-', ...
    'MarkerSize', 8, 'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', [1 0.7 0.7], 'LineWidth', 2);

% Plot Adults
errorbar(Week(1:length(Ldata)), Adata, CI.Adult, 'ko-', ...
    'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', [0.5 0.5 0.5], 'Color', 'k', 'LineWidth', 2);

% Axes and Labels
xlim([Week(1), Week(end)]);
ylabel("Tribolium Beetle Data", 'FontName', 'Times New Roman');
xlabel("Week", 'FontName', 'Times New Roman');

% Legend
lg = legend("Lurvae", "Pupae", "Adult", 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman');

hold off
