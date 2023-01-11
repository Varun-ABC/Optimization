function  [stress,sd] = standard_dev(n, r_in, r_out)
% n is number of quadruture points

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
elseif n == 5
    point = [-2.02018287, -0.9585724646, 0, 0.9585724646, 2.02018287];
    w = [0.01995324206, 0.3936193232,0.9453087205,0.3936193232, 0.01995324206]/ sqrt(pi);
elseif n == 6
    point = [-2.350604974, -1.335849074, -0.4360774119, 0.4360774119, 1.335849074, 2.350604974];
    w = [0.00453000991, 0.1570673203, 0.7246295952, 0.7246295952, 0.1570673203, 0.00453000991]/ sqrt(pi);
end


Nelem = 35;
L = 7.5; % semi-span in meters
E = 70e9; % Young's modulus, Pa
W = 0.5*500*9.8; % half of the operational weight, N
x = [0:L/Nelem:L].';
force = (2*(2.5*W)/(L^2))*[L:-L/Nelem:0].'; % nominal loading
% precomputing the cosine functions
c1 = cos((2-1)* pi * x / (2*L));
c2 = cos((2*2-1)* pi * x / (2*L));
c3 = cos((2*3-1)* pi * x / (2*L));
c4 = cos((2*4-1)* pi * x / (2*L));
c = [c1,c2,c3,c4];

% compute the displacements and the stresses

Iyy = CalcSecondMomentAnnulus(r_in, r_out);
stress = zeros(Nelem+1,1);
stress_sq = zeros(Nelem+1,1);
f0 = force(1);
% precomputing the stnadard devations
sig1 =  f0 / (10 * 1);
sig2 =  f0 / (10 * 2);
sig3 =  f0 / (10 * 3);
sig4 =  f0 / (10 * 4);
for i1 = 1:n
    p1 = point(i1) *(2)^.5 * sig1; % scale the points with the standared devation
    % since the mean is 0 it is not added here. 
    for i2 = 1:n
        p2 = point(i2) *(2)^.5 * sig2;
        for i3 = 1:n
            p3 = point(i3) *(2)^.5 * sig3;
            for i4 = 1:n
                p4 = point(i4) *(2)^.5 * sig4;
                delta_f = p1* c(:,1) + p2 *c(:,2) + p3 * c(:,3) + p4*c(:,4); 
                force_p = force + delta_f; % find the perturbed force
                u = CalcBeamDisplacement(L, E, Iyy, force_p, Nelem);
                s = CalcBeamStress(L, E, r_out, u, Nelem);
                stress = stress + w(i1) * w(i2) * w(i3)* w(i4) * s; % quadriture for the mean 
                stress_sq = stress_sq + w(i1) * w(i2) * w(i3)* w(i4) * s.^2;% quadriture for the mean squared
                % used in finding the standard devation
            end
        end
    end
end
sd = (stress_sq - stress.^2).^.5;
end
