function [cineq, ceq, jineq, jec] = noncon(RT, L, E, force, Nelem, Uts)
% finds non linear constrints and their jacobians
% Inputs:
% RT: Column vector of inner radii followed by thickness at each node
% L: Length of spar
% E: elastic modulus of material
% Downwards force at each node
% Nelem: number of elements 
% Uts: Ultimate tensile strength
% cineq: non linear inequality constriant
% ceq: non linear equality constraint 
% jineq: jacobian for the non linear ineq constraint
% jec: jacobian for the non linear eq constraint

    Nnode = Nelem + 1;    
    [cineq, ceq] = subcon(RT); %runs subfunction
    jineq = zeros(2* Nnode, Nnode);
    h = 1e-30; % complex step, step size
    for j = 1: 2 * Nnode
        % fills the jineq matrix
        % uses complex step to find the gradient of the stress at each
        % node with respect to the inner radius and the thickness
        xc = RT;
        xc(j) = xc(j) + complex(0,h);
        jineq(j,:) = imag(subcon(xc))/h;
    end
    jec = []; % no equality contraint
    function [cineq, ceq] = subcon(R_T)
        % finds non linear constrints
        % R_T is RT from above
        % outputs:
        % cineq, and ceq are same as above
        Iyy = zeros(Nelem + 1, 1);
        zmax = zeros(Nelem + 1, 1);
        for i = 1:Nnode
            % Calculates the moment of inetia of the spar at each node 
            % and the puts the outter radius in its own matrix
            ri = R_T(i);
            ro = ri + R_T(Nnode + i);
            zmax(i) = ro;
            Iyy(i) = pi / 4 *(ro^4 - ri^4);
            if real(Iyy(i)) < 0.0
                Iyy(i) = -Iyy(i);
            end
        end
        u = CalcBeamDisplacement(L, E, Iyy, force, Nelem);
        sigma = CalcBeamStress(L, E, zmax, u, Nelem);
        cineq = (sigma / Uts) - 1; % Ensusres that the stress will stay 
        % below the ultimate tensile strength
        ceq = []; % no equality constraint 
    end
end