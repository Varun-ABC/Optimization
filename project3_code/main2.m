% runs the hybrid optimization model 
clearvars;
close all;
% get the GPML Library

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
% read the sample data

samp = 1000;
w_lhs = xlsread('latin_samples_p3.xlsx',"A1:A1000");
r2_lhs = xlsread('latin_samples_p3.xlsx',"B1:B1000");
a1_lhs = xlsread('latin_samples_p3.xlsx',"C1:C1000");
st_dev = xlsread('latin_samples_p3.xlsx',"D1:D1000");
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
hyp1 = hyp;
% create two surrogates, one for each method
[f] = @(z) 1/ gp(hyp, @infExact, [], covfunc, likfunc, x, st_dev, z);
[f1] = @(z) 1/ gp(hyp1, @infExact, [], covfunc, likfunc, x, st_dev, z);
% using active- set
options = optimset("algorithm", "active-set"); %,'Display', 'iter');%, "PlotFcn","optimplotfirstorderopt" );
% define bounds
w_l = 3 * 2 * pi / 60; r2_l = .1; a1_l= 0;
w_u = 8 * 2 * pi / 60; r2_u = 1.5; a1_u= .3;
ub = [w_u, r2_u, a1_u];
lb = [w_l, r2_l, a1_l];
s = samp + 1;
max_sd = 0;
max_sd1 = 0;
%random initital condition
z1 = [(rand(1)+ 3 * 2 *pi/60)*5*2*pi/60, (rand(1) + .1) * 1.4, rand(1)*.3];
x1 = x;
st_dev1= st_dev;
for m = 1:10
z0 = [(rand(1)+ 3 * 2 *pi/60)*5*2*pi/60, (rand(1) + .1) * 1.4, rand(1)*.3];
% x = fmincon      (fun,x0,A,b, Aeq,beq,lb,ub, nonlcon,options)
[opt, sd] = fmincon(f, z0,[],[], [], [], lb, ub, [], options);
[opt1, sd1] = fmincon(f1, z1,[],[], [], [], lb, ub, [], options);
% opt = [w, r2, a1]
Tau = 3 * opt(2) * T;
Tspin = .1 * Tau;
fun = @(t,y) ode_(opt(1), opt(2), opt(3), y, t, r1);
[t,y_2] = ode45(fun, [Tspin,Tau],[0,0]);

fun1 = @(t,y) ode_(opt1(1), opt1(2), opt1(3), y, t, r1);
[t1,y_21] = ode45(fun1, [Tspin,Tau],[0,0]);

st_ = st_d(t,y_2, Tau, opt(1));
st_1 = st_d(t1,y_21, Tau, opt1(1));

if st_ > max_sd
    % check if theres a new max
    max_sd = st_;
    
end

if st_1 > max_sd1
    % check if theres a new max
    max_sd1 = st_1;
end
x(s,:) =  [opt(1), opt(2), opt(3)];
st_dev(s) = st_;
hyp = minimize(hyp, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
[f] = @(z) 1/ gp(hyp, @infExact, [], covfunc, likfunc, x, st_dev, z);
%
x1(s,:) =  [opt1(1), opt1(2), opt1(3)];
st_dev1(s) = st_1;
hyp1 = minimize(hyp1, @gp, -100, @infExact, [], covfunc,...
    likfunc, x1, st_dev1);
[f1] = @(z) 1/ gp(hyp1, @infExact, [], covfunc, likfunc, x1, st_dev1, z);
if st_ > st_1
    % if the random guess is greater than the gradual increase, set that to
    % new inital guess
    z1 = opt;
else 
    %otherwise keep recirculating
    z1 = opt1;
end
s = s+1;
end

