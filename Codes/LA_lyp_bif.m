%--------------------------------------------------------------------------
% 08/07/24, K.K., initial version
% 04/16/25, K.K., update version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bifurcation and Lyapunov exponent diagrams for LPAA model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc
tic

%% Global Plot Settings
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaulttextInterpreter', 'latex', 'defaultAxesFontSize', 16); 
format long

% Custom color definitions
bl = '#0072BD';     % Blue
yl = '#EDB120';     % Yellow
or = '#D95319';     % Orange

%% Simulation Parameters
back_step = 50;     % How many steps from end to record in bifurcation
no_iter = 1000;     % Total number of iterations

% Parameters from fitting (excluding first 4 time points)
b     = 39.4293;
c1    = 0.0317;
c2    = 0.0004;
mu_a  = 0.1287;
mu_l  = 0.4348;
A0    = 20;
L0    = 0;

% Pack parameters: [b, c1, c2, mu_a, mu_l, L0, A0]
best_params = [b, c1, c2, mu_a, mu_l, L0, A0];

%% Bifurcation Diagram Setup

bn = 100;                                      % Number of parameter steps
parameter_vec = linspace(0, 1, bn);            % Range of bifurcation parameter mu_a
L_plot = zeros(bn, back_step + 1);             % Store late L values
A_plot = zeros(bn, back_step + 1);             % Store late A values

disp("Creating bifurcation diagrams...")
for i = 1:length(parameter_vec)
    % Vary mu_a
    b_temp = parameter_vec(i);
    best_params(4) = b_temp;                   % Change mu_a in parameter vector

    % Simulate system
    [L, A, Jt] = LA(best_params, no_iter);

    % Save last `back_step + 1` values for plotting
    L_plot(i, :) = L(end - back_step:end);
    A_plot(i, :) = A(end - back_step:end);    
end

%% Lyapunov Exponent Computation

lyap = zeros(length(parameter_vec), 1);

disp("Generating Lyapunov diagrams...")
for j = 1:length(parameter_vec)
    disp("j = " + j)
    best_params(4) = parameter_vec(j);         % Vary mu_a again

    for k = 1:(no_iter - 1)
        % Iterate map (discard trajectory)
        [L, A, Jt] = LA(best_params, no_iter);
    end

    % Lyapunov exponent from Jacobian product norm
    lyap(j) = (1 / no_iter) * log(norm(Jt, "fro"));
end

toc

%% Plot Results

n = length(L) - 50:1:length(L);                % Not directly used in plot

figure()
t = tiledlayout(1, 2);

% Bifurcation Plot
nexttile
plot(parameter_vec, L_plot + A_plot, 'bl-');
ylabel('Total final Population')
xlabel('$\mu_a$')
pbaspect([2 1 1])
ax = gca;
ax.LineWidth = 1.5;
ax.FontWeight = 'bold';
ax.FontSize = 12;

% Lyapunov Exponent Plot
nexttile
plot(parameter_vec, lyap, 'bl.-', ...
     parameter_vec, zeros(size(parameter_vec)), 'r-', "LineWidth", 1.5);
ylabel("Lyapunov exponent")
xlabel("$\mu_a$")
pbaspect([2 1 1])
ax = gca;
ax.LineWidth = 1.5;
ax.FontWeight = 'bold';
ax.FontSize = 12;

% Optional Title
% title(t, "Arbitrary parameterization", "interpreter", "latex", "fontsize", 18)

%% Model Function Definitions

%--------------------------------------------------
% LA: LPAA Model Simulation
%--------------------------------------------------
function [L, A, Jlyap] = LA(params, Nmax)

    b     = params(1);
    c1    = params(2);
    c2    = params(3);
    mu_a  = params(4);
    mu_l  = params(5);

    L = zeros(Nmax, 1);
    A = zeros(Nmax, 1);

    L(1) = params(6);
    A(1) = params(7);

    Jlyap = evalJ(params, L(1), A(1));          % Initialize Jacobian

    for n = 1:Nmax-1
        % Update equations (discrete map)
        L(n+1) = A(n)*b*exp(-c1*A(n)) + (1 - mu_l)*L(n);  
        % L(n+1) = A(n)*b*exp(-c1*A(n) - c2*L(n)) + (1 - mu_l)*L(n); 
        A(n+1) = mu_l * L(n) * exp(-c2*A(n)) + A(n) * (1 - mu_a);

        % Accumulate Jacobian product for Lyapunov exponent
        J_temp = evalJ(params, L(n+1), A(n+1));
        Jlyap = J_temp * Jlyap;
    end
end

%--------------------------------------------------
% evalJ: Jacobian Evaluation at (L,A)
%--------------------------------------------------
function matrix = evalJ(params, L, A)

    b     = params(1);
    c1    = params(2);
    c2    = params(3);
    mu_a  = params(4);
    mu_l  = params(5);

    g1 = exp(-c1*A);               % Derivative w.r.t A for L(n+1)
    g2 = exp(-c2*A);               % Derivative w.r.t A for A(n+1)

    matrix = [1 - mu_l,     b * (1 - c1*A) * g1;
              mu_l * g2,   -c2 * mu_l * L * g2 + (1 - mu_a)];
end
