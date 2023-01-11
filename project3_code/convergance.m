% samples convergance 
clearvars
close all

mydir = 'gpml-matlab-v3.6-2015-07-07/';
addpath(mydir(1:end-1))
addpath([mydir,'cov'])
addpath([mydir,'doc'])
addpath([mydir,'inf'])
addpath([mydir,'lik'])
addpath([mydir,'mean'])
addpath([mydir,'prior'])
addpath([mydir,'util'])

w_lhs = xlsread('latin_samples_p3.xlsx',"A1:A2850");
r2_lhs = xlsread('latin_samples_p3.xlsx',"B1:B2850");
a1_lhs = xlsread('latin_samples_p3.xlsx',"C1:C2850");
st_dev = xlsread('latin_samples_p3.xlsx',"D1:D2850");
x = zeros(size(a1_lhs,1),3);
x(:,1) = w_lhs; x(:,2) = r2_lhs; x(:,3) = a1_lhs; 
r2 = 0.8;
a1 = 0.058;
r1 = 4.3;
T = 10000;
% noise convergance
% set the squared exponential covariance function
covfunc =  {@covMaterniso, 1}; % @covSEiso; 
% set the likelihood function to Gaussian
likfunc = @likGauss;
sn = 5; %1e-16; % this is the noise level

hyp_1.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)
hyp_2.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)
hyp_3.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)
hyp_4.cov = [log(1); log(.5)]; % first component is log(l) and second is log(sigma)

hyp_1.lik = log(sn);
hyp_2.lik = log(sn);
hyp_3.lik = log(sn);
hyp_4.lik = log(sn);


% maximize the likelihood function to find the hyperparameters
hyp_1 = minimize(hyp_1, @gp, -100, @infExact, [], covfunc,...
    likfunc, x(1:100,:), st_dev(1:100));
hyp_2 = minimize(hyp_2, @gp, -100, @infExact, [], covfunc,...
    likfunc, x(1:500,:), st_dev(1:500));
hyp_3 = minimize(hyp_3, @gp, -100, @infExact, [], covfunc,...
    likfunc, x(1:1000,:), st_dev(1:1000));
hyp_4 = minimize(hyp_4, @gp, -100, @infExact, [], covfunc,...
    likfunc, x(1:2850,:), st_dev(1:2850));
% create a diffrent surrogate for each, extract only the samples apprprate 
[f_100] = @(z) gp(hyp_1, @infExact, [], covfunc, likfunc, x(1:100,:), st_dev(1:100), z);
[f_500] = @(z) gp(hyp_2, @infExact, [], covfunc, likfunc, x(1:500,:), st_dev(1:500), z);
[f_1000] = @(z) gp(hyp_3, @infExact, [], covfunc, likfunc, x(1:1000,:), st_dev(1:1000), z);
[f_2850] = @(z) gp(hyp_4, @infExact, [], covfunc, likfunc, x(1:2850,:), st_dev(1:2850), z);

iter = .01;
sd_ode = zeros(length(.3:iter:.9),1);
sd_surr_1 = zeros(length(.3:iter:.9),1);
sd_surr_2 = zeros(length(.3:iter:.9),1);
sd_surr_3 = zeros(length(.3:iter:.9),1);
sd_surr_4 = zeros(length(.3:iter:.9),1);

i =1;
for w = .3:iter:.9
    % perform a sweep of the angular velocities
    fun = @(t,y) ode_(w, r2, a1, y, t, r1);
    Tau = 3 * w * T;
    [t,y_2] = ode45(fun, [0,Tau],[0,0]);
    y = y_2(:,1);
    y_m = 1/Tau * trapz(t,y);
    sd_ode(i) = 3* w * (1/Tau * trapz(t,(y - y_m).^2)).^(1/2);
    sd_surr_1(i) = f_100([w, r2, a1]);
    sd_surr_2(i) = f_500([w, r2, a1]);
    sd_surr_3(i) = f_1000([w, r2, a1]);
    sd_surr_4(i) = f_2850([w, r2, a1]);
    i = i+1;
end
% plot
plot(.3:iter:.9, sd_ode, 'DisplayName', 'ODE')
hold on
plot(.3:iter:.9, sd_surr_1, 'DisplayName', '100 Samples')
hold on
plot(.3:iter:.9, sd_surr_2, 'DisplayName', '500 Samples')
hold on
plot(.3:iter:.9, sd_surr_3, 'DisplayName', '1000 Samples')
hold on
plot(.3:iter:.9, sd_surr_4, 'DisplayName', '2850 Samples')
legend
xlabel('Angular Velocity')
ylabel('Standard Devation')