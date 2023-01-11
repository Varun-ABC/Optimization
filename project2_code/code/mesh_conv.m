clearvars
close all 
% This script performs a convergance study
% all the parameters are defined here
L = 7.5; %meter
rho = 1600; %kg/ m^3
E = 70e9; %Pa
Uts = 600e6; %Pa
mass = 500; %kg, max mass of the UAV
weight = 9.81 * mass; %N
Nelem = 30;
rmin = .01;
rmax = .05;
tmin = .0025;
elems = 5:5:70;
% the convergance code is run here:
masses = zeros(size(elems,2));
sigma_end = zeros(size(elems,2));
fmax = weight * 2.5 / L;

count = 0;
for i = elems
    count = count + 1;
    Nnode = i + 1;
    [RT,mass] = main(L, rho, E, Uts, weight, i, rmin, rmax, tmin);
    masses(count) = mass;
    force = linspace(fmax, 0, Nnode)';
    zmax = zeros(Nnode,1);
    Iyy = zeros(Nnode, 1);
    for j = 1:Nnode
        % generates a matrix of moment of inertia
        ri = RT(j);
        ro = ri + RT(Nnode + j);
        zmax(j) = ro;
        Iyy(j) = pi / 4 *(ro^4 - ri^4);
    end
    [u] = CalcBeamDisplacement(L, E, Iyy, force, i);
    [sigma] = CalcBeamStress(L, E, zmax, u, i);
    sigma_end(count) = sigma(Nnode)/sigma(1);
end

figure(1);
plot(elems, masses)
xlabel('Number of Elements')
ylabel('Mass of spar (kg)')


figure(2)
plot(elems,abs(sigma_end))
xlabel('Number of Elements')
ylabel('Normal stress at tip normalized to the root')

