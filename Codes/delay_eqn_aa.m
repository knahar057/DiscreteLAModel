%--------------------------------------------------------------------------
% 04/14/24, K.K., initial version
% 04/16/25, K.K., update version

% Parameter initialization
Nmax = 80;             % Maximum number of time steps
b = 52;                % Birth rate
c1 = 0.054;            % Nonlinear density-dependent effect (1)
c2 = 3*10^-4;          % Nonlinear density-dependent effect (2)
mu_a = 0.2;            % Adult mortality rate
mu_l = 0.3;            % Larval maturation rate

% Preallocate population vectors
L = zeros(Nmax, 1);    % Larval population (not updated in this code)
A = zeros(Nmax, 1);    % Adult population

% Initial values for adult population A(0)
aa = [10, 20, 40, 100, 200, 400];  % Different initial conditions

% Set initial larval population
L(1) = 0;

% Equilibrium value for adult population (analytical)
equib_A = log(b / mu_a) / (c1 + c2);  

% Initialize legend strings for plots
LegendStrings = cell(length(aa), 1);  

% Loop over different initial conditions
for k = 1:length(aa)
    A(1) = aa(k);  % Set initial adult population
    % Compute A(2) using given model dynamics
    A(2) = mu_l * L(1) * exp(-c2 * A(1)) + (1 - mu_a) * A(1);  
    
    % Iteratively compute A(n+1) using the delay equation model
    for n = 2:(Nmax - 1)
        A(n+1) = (1 - mu_a) * A(n) + ...
                 mu_l * b * A(n-1) * exp(-c1 * A(n-1) - c2 * A(n)) + ...
                 (1 - mu_l) * (A(n) - (1 - mu_a) * A(n-1)) * exp(c2 * (A(n-1) - A(n)));
    end
    
    % Plot each trajectory
    hold on
    plot(1:Nmax, A(1:Nmax), '-', 'LineWidth', 2)
    
    % Save label for legend
    LegendStrings{k} = ['A(0) = ', num2str(aa(k))];
end

% Plot the equilibrium line
plot(1:Nmax, equib_A * ones(1, Nmax), 'm*-')

% Add plot titles and labels
title('Adult population for delay equation')
ylabel('Adult', 'FontName', 'Times New Roman')
xlabel('Time', 'FontName', 'Times New Roman')

% Add legend
legend(LegendStrings, 'Interpreter', 'none')

hold off
