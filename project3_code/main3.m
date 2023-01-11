% runs brute force method
clearvars;
close all;

mydir = 'gpml-matlab-v3.6-2015-07-07/';
addpath(mydir(1:end-1))
addpath([mydir,'cov'])
addpath([mydir,'doc'])
addpath([mydir,'inf'])
addpath([mydir,'lik'])
addpath([mydir,'mean'])
addpath([mydir,'prior'])
addpath([mydir,'util'])
r1 = 4.3;
% dx = ode_(w, r2, a1, y, t, rho, m);
samp = 800;
w_lhs = xlsread('latin_samples_p3.xlsx',"A1:A800");
r2_lhs = xlsread('latin_samples_p3.xlsx',"B1:B800");
a1_lhs = xlsread('latin_samples_p3.xlsx',"C1:C800");
st_dev = xlsread('latin_samples_p3.xlsx',"D1:D800");
T = 10000;
x(:,1) = w_lhs; x(:,2) = r2_lhs; x(:,3) = a1_lhs; 
% set the squared exponential covariance function
covfunc = {@covMaterniso, 1}; 
hyp.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)

% set the likelihood function to Gaussian
likfunc = @likGauss;
sn = .005; %1e-16; % this is the noise level
hyp.lik = log(sn);

% maximize the likelihood function to find the hyperparameters
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
[f] = @(z) 1/ gp(hyp, @infExact, [], covfunc, likfunc, x, st_dev, z);

options = optimset("algorithm", "active-set"); %,'Display', 'iter');%, "PlotFcn","optimplotfirstorderopt" );
w_l = 3 * 2 * pi / 60; r2_l = .1; a1_l= 0;
w_u = 8 * 2 * pi / 60; r2_u = 1.5; a1_u= .3;
ub = [w_u, r2_u, a1_u];
lb = [w_l, r2_l, a1_l];
s = samp + 1;
max_sd = 0;
for m = 1:30
z0 = [(rand(1)+ 3 * 2 *pi/60)*5*2*pi/60, (rand(1) + .1) * 1.4, rand(1)*.3];
% x = fmincon      (fun,x0,A,b, Aeq,beq,lb,ub, nonlcon,options)
[opt, sd] = fmincon(f, z0,[],[], [], [], lb, ub, [], options);
% opt = [w, r2, a1]
Tau = 3 * opt(2) * T;
Tspin = .1 * Tau;
fun = @(t,y) ode_(opt(1), opt(2), opt(3), y, t, r1);
[t,y_2] = ode45(fun, [Tspin,Tau],[0,0]);

st_ = st_d(t,y_2, Tau, opt(1));

if st_ > max_sd
    max_sd = st_;
end

x(s,:) =  [opt(1), opt(2), opt(3)];
st_dev(s) = st_;
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
[f] = @(z) 1/ gp(hyp, @infExact, [], covfunc, likfunc, x, st_dev, z);
s = s+1;
end

