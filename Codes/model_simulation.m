%--------------------------------------------------------------------------
% 03/24/24, K.K., initial version
% 04/16/25, K.K., update version

% Parameter initialization
Nmax = 100;               % Maximum number of time steps
b = 52;                   % Birth rate (case: b > mu_a)
c1 = 0.054;               % Density-dependent factor (1)
c2 = 3 * 10^-4;           % Density-dependent factor (2)
mu_a = 0.2;               % Adult mortality rate
mu_l = 0.3;               % Larval maturation rate

% Preallocate population vectors
L = zeros(Nmax, 1);       % Larval population
A = zeros(Nmax, 1);       % Adult population

% Initial conditions
L(1) = 0;                 % Initial larval population
A(1) = 20;                % Initial adult population

%% Compute equilibrium points (only valid when b > mu_a)
equib_A = log(b / mu_a) / (c1 + c2);                             % Adult equilibrium
equib_L = (b / mu_l) * equib_A * exp(-c1 * equib_A);            % Larval equilibrium

%% Main time-stepping loop
for n = 1:(Nmax - 1)
    % Update larval population (density-dependent birth + survival)
    L(n+1) = A(n) * b * exp(-c1 * A(n)) + (1 - mu_l) * L(n); 
    
    % Alternate form if L feedback included:
    % L(n+1) = A(n) * b * exp(-c1 * A(n) - c2 * L(n)) + (1 - mu_l) * L(n);

    % Update adult population (from larvae maturation + adult survival)
    A(n+1) = mu_l * L(n) * exp(-c2 * A(n)) + (1 - mu_a) * A(n);
end

% Plotting results
figure()
title('b > \mu_a')  % Modify to 'b < \mu_a' if testing that regime
hold on

plot(1:Nmax, A(1:Nmax), 'g.-', 'LineWidth', 2)     % Adult population
plot(1:Nmax, L(1:Nmax), 'b.-', 'LineWidth', 2)     % Larval population

% Plot equilibrium levels as horizontal lines
plot(1:Nmax, equib_A * ones(1, Nmax), 'r*-')       
plot(1:Nmax, equib_L * ones(1, Nmax), 'r*-')       

% Axis labels
ylabel("Adult & Larvae", 'FontName', 'Times New Roman')
xlabel("Time", 'FontName', 'Times New Roman')

% Legend for b > mu_a case
legend("A", "L", "Equilibrium", 'FontName', 'Times New Roman')

hold off
