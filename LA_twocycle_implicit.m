%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes a figure showing the existence of a two-cycle using the
% discrete delay form of the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

% setup tiledlayout
t = tiledlayout(2,3,Padding="tight");

fs = 24;


%% implicit plot c2 = 0
% needs to be in the form f(x,y) = 0 for plotting.

b = 19;
mu_l = 0.54; s_l = (1-mu_l);
mu_a = 0.95; s_a = (1-mu_a);
c1 = 0.18;
c2 = 0;
param = [b mu_l mu_a c1 c2];

% functions
f1 = @(p,q) -p + b*(1-s_l)*p.*exp(-c1*p - c2*q)...
        + s_l.*exp(c2*(p-q)).*(q - s_a*p) + s_a*q;

f2 = @(p,q) -q + b*(1-s_l)*q.*exp(-c1*q - c2*p)...
        + s_l.*exp(c2*(q-p)).*(p - s_a*q) + s_a*p;

nexttile(1)
    fimplicit(f1,[-10,200,-10,200],'r', 'LineWidth',2); hold on
    fimplicit(f2,'b','LineWidth',2);
    fimplicit(@(p,q)p-q, 'k--','LineWidth',1)
    set(gca,'FontSize',12)

    xlabel('$p$','Interpreter','latex','FontSize',fs);
    ylabel('$q$','Interpreter','latex','FontSize',fs)

    xlim([0 30]); ylim([0 30])
    
    axis square
    title('(a) $c_2$ = 0','Interpreter','latex','FontSize',fs)
    legend({'Eqn. 19','Eqn. 20','$p=q$'},'interpreter','latex', 'Location','northwest','FontSize',fs)

    ax=gca;
    ax.FontSize = 20

A1 = run_A(param,500);

nexttile(4)
    plot(A1,'k','LineWidth',1,'Marker','.','MarkerSize',15);
    xlim([450 500]); ylim([0 25])
    set(gca,'FontSize',12)
    xlabel(t, 'Iteration $t$', 'Interpreter','latex','FontSize',fs)
    ylabel('Adults', 'Interpreter','latex','FontSize',fs)
    axis square
    ax=gca;
    ax.FontSize = 20

%% implicit plot c1 = 0.163
b = 19;
mu_l = 0.54; s_l = (1-mu_l);
mu_a = 0.95; s_a = (1-mu_a);
c1 = 0.18;
c2 = 0.163;

param = [b mu_l mu_a c1 c2];

% functions
f1 = @(p,q) -p + b*(1-s_l)*p.*exp(-c1*p - c2*q)...
        + s_l.*exp(c2*(p-q)).*(q - s_a*p) + s_a*q;

f2 = @(p,q) -q + b*(1-s_l)*q.*exp(-c1*q - c2*p)...
        + s_l.*exp(c2*(q-p)).*(p - s_a*q) + s_a*p;

nexttile(2)
    fimplicit(f1,[-10,200,-10,200],'r', 'LineWidth',2); hold on
    fimplicit(f2,'b','LineWidth',2);
    fimplicit(@(p,q)p-q, 'k--','LineWidth',1)
    set(gca,'FontSize',12)

    xlabel('$p$','Interpreter','latex','FontSize',fs);
    ylabel('$q$','Interpreter','latex','FontSize',fs)

    xlim([0 30]); ylim([0 30])
    
    axis square
    title('(b) $c_2$ = 0.15','Interpreter','latex','FontSize',fs)
    ax=gca;
    ax.FontSize = 20

A2 = run_A(param,500);

nexttile(5)
    plot(A2,'k','LineWidth',1,'Marker','.','MarkerSize',15);
    xlim([450 500]); ylim([0 25])
    set(gca,'FontSize',12)
    axis square
    ax=gca;
    ax.FontSize = 20

%% implicit plot c2 = 0.30
b = 19;
mu_l = 0.54; s_l = (1-mu_l);
mu_a = 0.95; s_a = (1-mu_a);
c1 = 0.18;
c2 = 0.30;

param = [b mu_l mu_a c1 c2];

% functions
f1 = @(p,q) -p + b*(1-s_l)*p.*exp(-c1*p - c2*q)...
        + s_l.*exp(c2*(p-q)).*(q - s_a*p) + s_a*q;

f2 = @(p,q) -q + b*(1-s_l)*q.*exp(-c1*q - c2*p)...
        + s_l.*exp(c2*(q-p)).*(p - s_a*q) + s_a*p;

nexttile(3)
    fimplicit(f1,[-10,200,-10,200],'r', 'LineWidth',2); hold on
    fimplicit(f2,'b','LineWidth',2);
    fimplicit(@(p,q)p-q, 'k--','LineWidth',1)
    set(gca,'FontSize',12)

    xlabel('$p$','Interpreter','latex','FontSize',fs);
    ylabel('$q$','Interpreter','latex','FontSize',fs)

    xlim([0 30]); ylim([0 30])
    
    axis square
    title('(c) $c_2$ = 0.30','Interpreter','latex','FontSize',fs)
    ax=gca;
    ax.FontSize = 20

A3 = run_A(param,500);

nexttile(6)
    plot(A3,'k-','LineWidth',1,'Marker','.','MarkerSize',15);
    xlim([450 500]);ylim([0 25])
    set(gca,'FontSize',12)
    axis square
    ax=gca;
    ax.FontSize = 20





function A = run_A(param,N)
    % parameters
    b = param(1);
    mu_l = param(2);
    s_l = (1-mu_l);
    mu_a = param(3);
    s_a = (1-mu_a);
    c1 = param(4);
    c2 = param(5);

    % initial conditions (trying to start near the cycle)
    A = NaN(1,N);
    A(1) = 4;
    A(2) = 10;
    
    % loop using delay difference equation
    for j = 2:N
        A(j+1) = b*mu_l*A(j-1)*exp(-c1*A(j-1) - c2*A(j))... 
                    + s_l*exp(c2*(A(j-1) - A(j)))*(A(j) - s_a*A(j-1))...
                    + s_a*A(j);
        
    end

end