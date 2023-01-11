function [a0, afinal, flux0, flux1, improval] = main()
    % takes no inputs all values are hard coded in
    % outputs: 
    % a0: inital guess for design variables
    % afinal: final matrix of design variables
    % flux0: initial flux
    % flux1: final flux
    % improval: percent change in flux
    close all
    L = .05;
    Nx = 100;
    Ny = 100;
    x = [0:L/Nx:L].';
    kappa = 20;
    Ttop = 20 + 273;
    Tbot= 90 + 273; % hard coded knowns
    fun = @(a) objection_func(a, L, Nx, Ny, kappa, Ttop, Tbot, x);
    %line above: an anonomus function that eliminates all other imputs
    %other than the one we want to optimize, the matrix of design variables
    a0 = [.03; 0.01.*rand(4,1)]; % inital guess, randomized
    a0 = [0.03; 0; 0; 0; 0.02; .0; 0; .00; 0];
    [Aineq, Bineq] = limits(a0, Nx, L, x); %find contstraint function
    options = optimset('Display', 'iter', 'algorithm', 'sqp',"PlotFcn", 'optimplotfirstorderopt');
    afinal = fmincon(fun,a0,Aineq,Bineq, [],[],[],[],[], options); 
    %line above: run fmincon from matlab optimization toolbox
    h = height(a0, L, x);% calculate the initial and final conditions 
    %and determine the improvement 
    h1 = height(afinal, L, x);
    flux0 = CalcFlux(L, h, Nx, Ny, kappa, Ttop, Tbot);
    flux1 = CalcFlux(L, h1, Nx, Ny, kappa, Ttop, Tbot);
    plot(x,h, x,h1), legend("Initial", "Final"), ylabel('Height (m)') , xlabel('Distance along X-axis (m)'), title('2D Cross section');
    improval = (flux1 - flux0) / flux0 * 100;