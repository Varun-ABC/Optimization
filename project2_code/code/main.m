function [opt, mass] = main(L, rho, E, Uts, weight, Nelem, rmin, rmax, tmin)
% MAIN function, defines the objective function and the bounderies and
% inequality constraints
% plugs them into fmincon optimization 
% Inputs: L: Length of the spar
% rho: density of the material
% L: Length of spar
% E: elastic modulus of material
% weight: weight of the whole UAV
% Nelem: number of elements 
% Uts: Ultimate tensile strength  
% rmin: minimum inner radius
% rmax: Max outter radius
% tmin: min thickness
% Outputs:
% opt: optimized radius and thickness
% Optimized mass

%function [w] = obj(RT, Nnode , L, rho)
%function [cineq, ceq, Jineq, Jeq] = noncon(RT, L, E, force, Nelem, Uts)
%function [Aineq, B, lb, ub] = bounds(Nnode, rmin, rmax, tmin)

    Nnode = Nelem + 1;
    fun = @(RT) obj(RT, Nnode , L, rho); % creates an objective function with only one input
    fmax = weight * 2.5 / L;
    force = linspace(fmax, 0, Nnode);%  the force distribution for this case is triangular
    % fmax is calculated using geometetry and we need to know thr force on each
    % node
    nonlin = @(RT) noncon(RT, L, E, force, Nelem, Uts);
    [Aineq, B, lb, ub] = bounds(Nnode, rmin, rmax, tmin);
    x0 = [linspace(rmax-.01,rmin,Nnode)'; ones(Nnode, 1)*.01];
    options = optimset("algorithm", "active-set",'Display', 'iter', ...
        "GradConstr", "on", "GradObj", "on", "DerivativeCheck", "off", "PlotFcn","optimplotfirstorderopt" );
    [opt, mass] = fmincon(fun, x0, Aineq, B, [],[], lb, ub, nonlin, options);
end