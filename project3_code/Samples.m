% this code generates samples and finds the value of the ODE at those
% sample locations

clearvars
close all 

samp = 2850;
% use Latin Hypercube to find samples
w_lhs = lhsdesign(samp,1);
% need to scale each set of samples to fit within its bounds and to get the
% full range of bounds
w_lhs = (w_lhs + 3 * 2 *pi/60)*5*2*pi/60;
r2_lhs = lhsdesign(samp,1);
r2_lhs = (r2_lhs + .1) * 1.4;
a1_lhs = lhsdesign(samp,1);
a1_lhs = a1_lhs * .3;
st_dev = zeros(samp,1);
T = 10000;

r1 = 4.3;

for i = 1:samp
    % solve the ode at each set of samples
    fun = @(t,y) ode_(w_lhs(i), r2_lhs(i), a1_lhs(i), y, t, r1);
    Tau = 3 * w_lhs(i) * T;
    Tspin = .1 * Tau;
    [t,y_2] = ode45(fun, [Tspin, Tau],[0,0]);
    y = y_2(:,1);
    y_m = 1/Tau * trapz(t,y);
    st_dev(i) = 3* w_lhs(i) * (1/Tau * trapz(t,(y - y_m).^2)).^(1/2);
end
x = zeros(samp,4);
x(:,1) = w_lhs; x(:,2) = r2_lhs; x(:,3) = a1_lhs; x(:,4) = st_dev;
% save to an excel file
xlswrite('latin_samples_p3.xlsx',x)
