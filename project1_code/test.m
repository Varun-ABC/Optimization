function out = test ()
% takes no inputs, all inputs are is hard coded in
% outputs the flux
% graphs the effect of increasing values for Nx and Ny
close all
a0 = [0.0300, -0.0016, -0.0199, 0.0017];
L = .05;
kappa = 20;
Ttop = 20 + 273;
Tbot= 90 + 273; % inputs to CalcFlux
start = 10; %starting Nx and Ny value
ending = 300; %ending Nx and Ny values
iter = 10;
n = ((ending - start)/ iter) +1; % find the total amount of-- 
%Nx and Ny values to be testes
x_flux = zeros(n,1);
y_flux = zeros(n,1);
Nx = 100;
Ny = 100;
counter  = 0;
x_axis = zeros(n,1);
for i = start:iter: ending
    % calculate flux for Nx and Ny values
    counter = counter + 1;
    x_axis(counter) = i;
    x_x = [0:L/i:L].';
    x_y = [0:L/Nx:L].';
    h_y = height(a0, L, x_y);
    h_x = height(a0, L, x_x);
    fx = CalcFlux(L, h_x, i, Ny, kappa, Ttop, Tbot);
    fy = CalcFlux(L, h_y, Nx, i, kappa, Ttop, Tbot);
    x_flux(counter) = fx(1); 
    y_flux(counter) = fy(1);
end
flux_max = zeros(n,1);
flux_max = flux_max + max(x_flux(n), y_flux(n));
plot(x_axis, x_flux,x_axis, y_flux, x_axis, flux_max * .99), 
xlabel('Nx and Ny'), ylabel('flux'), 
legend('Flux changing Nx','Flux changing Ny', '99% of maximum value'), title('Mesh Convergance')
out = flux_max;