function [f] = objection_func(a, L, Nx, Ny, kappa, Ttop, Tbot, x)
% this wrapper function calculates the flux using all the givens in the
% problem statement as well as the matrix of design variables
    h = height(a, L, x); %need h to calculate the flux
    [flux,~,~,~] = CalcFlux(L, h, Nx, Ny, kappa, Ttop, Tbot);%supress all 
    %other outputs other than flux
    f = 1/flux; %b/c fmincon is going to find the min we need to give it 
    % the reciprocal 