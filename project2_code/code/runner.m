clearvars
close all 
% This script runs the optimization code
% all the parameters are defined here
L = 7.5; %meter
rho = 1600; %kg/ m^3
E = 70e9; %Pa
Uts = 600e6; %Pa
mass = 500; %kg, max mass of the UAV
weight = 9.81 * mass; %N
Nelem = 35;
rmin = .01;
rmax = .05;
tmin = .0025;
% the optimization code is run here:
[RT,mass] = main(L, rho, E, Uts, weight, Nelem, rmin, rmax, tmin);
% RT is the optimized inner radius and thicknesses at each node; mass is the optimized mass
ri = RT(1:Nelem+1); % seperates the column matrix for the inner radius
ro = RT(1:Nelem+1) + RT(Nelem+2:2*(Nelem+1)); % generates the column matrix 
% for the outer matrix by adding the thickness to the inner radius
x = linspace(0,L, Nelem+1);
Nnode = Nelem + 1;
figure(2);
plot(x, ri)
hold on
plot(x,ro)
xlabel('distance along wing (m)')
ylabel('Radii (m)')
legend("inner radius","outter radius")
fmax = weight * 2.5 / L;
force = linspace(fmax, 0, Nnode)'; % the force distribution for this case is triangular
% fmax is calculated using geometetry and we need to know thr force on each
% node
zmax = zeros(Nnode,1);
Iyy = zeros(Nelem + 1, 1);
weight = mass  * 9.81;
    for i = 1:Nnode
        % generates a matrix of moment of inertia
        ri = RT(i);
        ro = ri + RT(Nnode + i);
        zmax(i) = ro;
        Iyy(i) = pi / 4 *(ro^4 - ri^4);
    end
[u] = CalcBeamDisplacement(L, E, Iyy, force, Nelem);
[sigma] = CalcBeamStress(L, E, zmax, u, Nelem);

figure(3)
plot(x,sigma,'ks-')
xlabel('distance along wing (m)')
ylabel('magnitude of normal stress (Pa)')

figure(4)
plot(x,u(2:2:2*Nnode),'ks-')
xlabel('distance along wing (m)')
ylabel('deflection of spar (m)')
% figure(4)
% plot(x,sigma/Uts,'ks-')
% xlabel('distance along wing (m)')
% ylabel('magnitude of normal stress normilized to root')

%function [w] = obj(RT, Nnode , L, rho)
%function [cineq, ceq, Jineq, Jeq] = noncon(RT, L, E, force, Nelem, Uts)
%function [Aineq, B, lb, ub] = bounds(Nnode, rmin, rmax, tmin)
