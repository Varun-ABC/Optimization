% This script runs convergance for noise

clearvars
close all
% get the GMPL library
mydir = 'gpml-matlab-v3.6-2015-07-07/';
addpath(mydir(1:end-1))
addpath([mydir,'cov'])
addpath([mydir,'doc'])
addpath([mydir,'inf'])
addpath([mydir,'lik'])
addpath([mydir,'mean'])
addpath([mydir,'prior'])
addpath([mydir,'util'])
% read samples from sheet
w_lhs = xlsread('latin_samples_p3.xlsx',"A1:A800");
r2_lhs = xlsread('latin_samples_p3.xlsx',"B1:B800");
a1_lhs = xlsread('latin_samples_p3.xlsx',"C1:C800");
st_dev = xlsread('latin_samples_p3.xlsx',"D1:D800");
sd_ode = xlsread('latin_samples_p3.xlsx',"E1:E61");
x = zeros(size(a1_lhs,1),3);
x(:,1) = w_lhs; x(:,2) = r2_lhs; x(:,3) = a1_lhs; 
r2 = 0.8;
a1 = 0.058;
r1 = 4.3;
T = 10000;
% set the squared exponential covariance function
covfunc =  {@covMaterniso, 1}; % @covSEiso; 
% set the likelihood function to Gaussian
likfunc = @likGauss;
hyp_1.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)
hyp_2.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)
hyp_3.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)
hyp_4.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)
hyp_5.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)

hyp_1.lik = log(1e-16); % this is the noise level
hyp_2.lik = log(.00005);
hyp_3.lik = log(.05);
hyp_4.lik = log(5);
hyp_5.lik = log(50);


% maximize the likelihood function to find the hyperparameters
hyp_1 = minimize(hyp_1, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
hyp_2 = minimize(hyp_2, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
hyp_3 = minimize(hyp_3, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
hyp_4 = minimize(hyp_4, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
hyp_5 = minimize(hyp_5, @gp, -100, @infExact, [], covfunc,...
    likfunc, x, st_dev);
% Find the surrogate models for each
[f_1] = @(z) gp(hyp_1, @infExact, [], covfunc, likfunc, x, st_dev, z);
[f_2] = @(z) gp(hyp_2, @infExact, [], covfunc, likfunc, x, st_dev, z);
[f_3] = @(z) gp(hyp_3, @infExact, [], covfunc, likfunc, x, st_dev, z);
[f_4] = @(z) gp(hyp_4, @infExact, [], covfunc, likfunc, x, st_dev, z);
[f_5] = @(z) gp(hyp_5, @infExact, [], covfunc, likfunc, x, st_dev, z);


iter = .01;
sd_surr_1 = zeros(length(.3:iter:.9),1);
sd_surr_2 = zeros(length(.3:iter:.9),1);
sd_surr_3 = zeros(length(.3:iter:.9),1);
sd_surr_4 = zeros(length(.3:iter:.9),1);
sd_surr_5 = zeros(length(.3:iter:.9),1);

i =1;
for w = .3:iter:.9
    sd_surr_1(i) = f_1([w, r2, a1]);
    sd_surr_2(i) = f_2([w, r2, a1]);
    sd_surr_3(i) = f_3([w, r2, a1]);
    sd_surr_4(i) = f_4([w, r2, a1]);
    sd_surr_5(i) = f_5([w, r2, a1]);
    i = i+1;
end
plot(.3:iter:.9, sd_ode, 'DisplayName', 'ODE')
hold on
plot(.3:iter:.9, sd_surr_1, 'DisplayName', '1e-16')
hold on
plot(.3:iter:.9, sd_surr_2, 'DisplayName', '.0005 noise level')
hold on
plot(.3:iter:.9, sd_surr_3, 'DisplayName', '.05 noise level')
hold on
plot(.3:iter:.9, sd_surr_4, 'DisplayName', '5 noise level')
hold on
plot(.3:iter:.9, sd_surr_5, 'DisplayName', '500 noise level')
legend
xlabel('Angular Velocity')
ylabel('Standard Devation')