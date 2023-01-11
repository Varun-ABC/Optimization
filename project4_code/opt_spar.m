% minimize wing spar weight subject to stress constraints at manuever
clearvars;
close all;

% carbon fiber values from http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
Nelem = 35;
L = 7.5; % semi-span in meters
rho = 1600; % density of standard carbon fiber, kg/m^3
yield = 600e6; % tensile strength of standard carbon fiber, Pa
E = 70e9; % Young's modulus, Pa
W = 0.5*500*9.8; % half of the operational weight, N
force_nom = (2*(2.5*W)/(L^2))*[L:-L/Nelem:0].'; % loading at manueuver
delta1 = 0;
x = [0:L/Nelem:L].';

% define function and constraints
fun = @(r) SparWeight(r, L, rho, Nelem);
% precomputinf cosine functions
c1 = cos((2-1)* pi * x / (2*L));
c2 = cos((2*2-1)* pi * x / (2*L));
c3 = cos((2*3-1)* pi * x / (2*L));
c4 = cos((2*4-1)* pi * x / (2*L));
c = [c1,c2,c3,c4];
n = 2;
nonlcon = @(r) WingConstraints(r, L, E, force_nom, yield, Nelem, c,n);
lb = 0.01*ones(2*(Nelem+1),1);
up = 0.05*ones(2*(Nelem+1),1);
A = zeros(Nelem+1,2*(Nelem+1));
b = -0.0025*ones(Nelem+1,1);
for k = 1:(Nelem+1)
    A(k,k) = 1.0;
    A(k,Nelem+1+k) = -1.0;
end

% define initial guess (the nominal spar)
r0 = zeros(2*(Nelem+1),1);
r0(1:Nelem+1) = 0.0415*ones(Nelem+1,1);
r0(Nelem+2:2*(Nelem+1)) = 0.05*ones(Nelem+1,1);

options = optimset('GradObj','on','GradConstr','on', 'TolCon', 1e-4, ...
    'TolX', 1e-8, 'Display','iter', 'algorithm', 'sqp', "PlotFcn","optimplotfirstorderopt"); %, 'DerivativeCheck','on');
[ropt,fval,exitflag,output] = fmincon(fun, r0, A, b, [], [], lb, up, ...
     nonlcon, options);

% plot optimal radii
r_in = ropt(1:Nelem+1);
r_out = ropt(Nelem+2:2*(Nelem+1));
figure
plot(x, r_in, '-ks');
hold on;
plot(x, r_out, '--ks');
ylabel('Radii (m)')
xlabel('distance along spar (m)')
% display weight and stress constraints
[f,~] = fun(ropt);
[c,~,~,~] = nonlcon(ropt);
Iyy = CalcSecondMomentAnnulus(r_in, r_out);
[u] = CalcBeamDisplacement(L, E, Iyy, force_nom, Nelem);
[sigma] = CalcBeamStress(L, E, r_out, u, Nelem);

figure
plot(x,sigma,'ks-')
xlabel('distance along wing (m)')
ylabel('magnitude of normal stress (Pa)')

figure
plot(x,u(1:2:2*Nelem+1),'ks-')
xlabel('distance along wing (m)')
ylabel('deflection of spar (m)')

figure 
[stress,sd] = standard_dev(2, r_in, r_out);
figure
plot(x, stress)
hold on
plot(x, stress + 6*sd);
hold on
plot(x, stress - 6*sd);
hold on
legend({'mean stress','mean + 6 sd','mean - 6 sd'},'Location', 'northeast')
xlabel('Distance along span (m)') 
ylabel('stress (Pa)')

