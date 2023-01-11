function dev = st_d(t,y_2, Tau, w)
% This function finds the standard devation of d_phi/d_tau
% inputs: 
% t: list of tau values
% y_2: list of  d_phi/d_tau and phi values
% Tau: total non-dimentional time to run for
% w: average angular velocity
% output:
% dev: standard devation
y = y_2(:,1);
y_m = 1/Tau * trapz(t,y); % this is the mean of the standard devation
dev = 3* w * (1/Tau * trapz(t,(y - y_m).^2)).^(1/2);