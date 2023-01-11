function [dx] = ode_(w, r2, a1, y, t, r1)
% is a representaion of the equasion found in the paper
% inputs:
% w: average angular velocity
% r2
% a1
% y = [y1, y2]
% y1 = phi 
% y2 = dphi / dtau
% t: simulation time
% r1
% outputs:
% dx = [d_phi/d_tau, d^2_phi/d_tau^2]

g = 9.81;
Q0 = 20;
a0 = 0.036;
tau = t;
gam = (1/(3*w)) * (g / r2)^.5;
alpha = a0 - a1* cos(tau);
beta = 3 * a1 * sin(tau);
eps = r1 / (9 * r2);
A = gam / Q0;
B = (eps - gam^2 * alpha);
C = gam^2 * beta;
% create a state vecor here
dx = [0,0]';
dx(1)= -A * y(1) - B * sin(y(2)) - C*cos(y(2));
dx(2)= y(1);
