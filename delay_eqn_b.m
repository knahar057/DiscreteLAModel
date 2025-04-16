%--------------------------------------------------------------------------
% 03/24/24, K.K., initial version
% 04/16/25, K.K., update version

% Parameter initialization
Nmax = 100;                          % Maximum number of time steps
b = [0.1, 0.2, 3, 10, 35, 52];       % Different birth rates
c1 = 0.054;                          % Nonlinear density-dependent effect (1)
c2 = 3 * 10^-4;                      % Nonlinear density-dependent effect (2)
mu_a = 0.2;                          % Adult mortality rate
mu_l = 0.3;                          % Larval maturation rate

% Preallocate population vectors
L = zeros(Nmax, 1);                  % Larval population (not updated in loop)
A = zeros(Nmax, 1);                  % Adult population

% Set initial conditions
L(1) = 0;                            % Initial larval population
A(1) = 20;                           % Initial adult population
A(2) = mu_l * L(1) * exp(-c2 * A(1)) + (1 - mu_a) * A(1);  % Compute A(2)

% Initialize legend strings
LegendsStrings = cell(length(b), 1);

% Loop over different birth rates
for k = 1:length(b)
    % Simulate dynamics for each b(k)
    for n = 2:(Nmax - 1)
        A(n+1) = (1 - mu_a) * A(n) + ...
                 mu_l * b(k) * A(n-1) * exp(-c1 * A(n-1) - c2 * A(n)) + ...
                 (1 - mu_l) * (A(n) - (1 - mu_a) * A(n-1)) * exp(c2 * (A(n-1) - A(n)));
    end
    
    % Plot the solution
    hold on
    plot(1:Nmax, A(1:Nmax), '-', 'LineWidth', 2)
    
    % Store legend label
    LegendsStrings{k} = ['b = ', num2str(b(k))];
end

% Plot equilibrium levels (only for b(3) and onward)
for k = 3:length(b)
    equib_A = log(b(k) / mu_a) / (c1 + c2);  % Analytical equilibrium
    plot(1:Nmax, equib_A * ones(1, Nmax), 'r.-')  % Horizontal line at equilibrium
end

% Add title and axis labels
title('Adult population for delay equation with \mu_a = 0.2')
ylabel('Adult', 'FontName', 'Times New Roman')
xlabel('Time', 'FontName', 'Times New Roman')

% Add legend
legend(LegendsStrings, 'Interpreter', 'none')
