
function [cineq, ceq, dcdx, dceqdx] = WingConstraints(r, L, E, force, yield, Nelem, c, n)
% Computes the nonlinear inequality constraints for the wing-spar problem
% Inputs:
%   r - the DVs; x(1:Nelem+1) inner and x(Nelem+2:2*(Nelem+1) outer radius
%   L - length of the beam
%   E - longitudinal elastic modulus
%   force - nominal force per unit length along the beam axis x
%   yield - the yield stress for the material
%   Nelem - number of finite elements to use
%   c: The Cosine functions
%   n: number of quad points to be used
% Outputs:
%   cineq, ceq - inequality (stress) and equality (empty) constraints
%   dcdx, dceqdx - Jacobians of cineq and ceq
%--------------------------------------------------------------------------
assert( size(force,1) == (Nelem+1) );
assert( size(r,1) == (2*(Nelem+1)) );

[cineq] = CalcInequality(r, n);
ceq = [];
dcdx = zeros(2*(Nelem+1),Nelem+1);
dceqdx = [];

for k = 1:2*(Nelem+1)
    xc = r;
    xc(k) = xc(k) + complex(0.0, 1e-30);
    xc = CalcInequality(xc, n);
    dcdx(k,:) = imag(xc)/1e-30;
end 

    function [cineq] = CalcInequality(r, n)
        if n == 1
            point = [0];
            w = [sqrt(pi)]/ sqrt(pi);
        elseif n == 2
            point = [-0.7071067811865475244008, 0.7071067811865475244008];
            w = [0.8862269254527580136491, 0.8862269254527580136491]/ sqrt(pi);
        elseif n==3
            point = [-1.224744871391589049099,0, 1.224744871391589049099];
            w = [0.295408975150919337883, 1.181635900603677351532, 0.295408975150919337883]/ sqrt(pi);
        elseif n == 4
            point = [-sqrt((3-sqrt(6))/2),sqrt((3-sqrt(6))/2), -sqrt((3+sqrt(6))/2), sqrt((3+sqrt(6))/2)];
            w = 1/sqrt(pi)*[sqrt(pi)/ (4 * (3- sqrt(6))), sqrt(pi)/ (4 * (3- sqrt(6))),...
            sqrt(pi)/ (4 * (3+ sqrt(6))), sqrt(pi)/ (4 * (3 + sqrt(6)))];
        end
        % compute the displacements and the stresses
        r_in = r(1:Nelem+1);
        r_out = r(Nelem+2:2*(Nelem+1));
        Iyy = CalcSecondMomentAnnulus(r_in, r_out);
        stress = zeros(Nelem+1,1);
        stress_sq = zeros(Nelem+1,1);
        f0 = force(1); % force at root
        sig1 =  f0 / (10 * 1); % the standard devations of each random variable
        sig2 =  f0 / (10 * 2);
        sig3 =  f0 / (10 * 3);
        sig4 =  f0 / (10 * 4);
        % four for loops for four random variables
        for i1 = 1:n
            p1 = point(i1) *(2)^.5 * sig1; % scale the points with the standared devation
    % since the mean is 0 it is not added here. 
            for i2 = 1:n
                p2 = point(i2) *(2)^.5 * sig2;
                for i3 = 1:n
                    p3 = point(i3) *(2)^.5 * sig3;
                    for i4 = 1:n
                        p4 = point(i4) *(2)^.5 * sig4;
                        delta_f = p1* c(:,1) + p2 *c(:,2) + p3 * c(:,3) + p4 * c(:,4);
                        force_p = force + delta_f; % find the perturbed force
                        u = CalcBeamDisplacement(L, E, Iyy, force_p, Nelem);
                        s = CalcBeamStress(L, E, r_out, u, Nelem);
                        stress = stress + w(i1) * w(i2) * w(i3)* w(i4) * s; % quadriture for the mean 
                        stress_sq = stress_sq + w(i1) * w(i2) * w(i3)* w(i4) * s.^2; % quadriture for the mean squared
                                    % used in finding the standard devation

                    end
                end
            end
        end
         sd = (abs(stress_sq - stress.^2)).^.5;
         cineq = (stress + 6 * sd)/yield - 1;
    end
end
