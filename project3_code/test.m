% This code is created to test the ODE to ensure it is working
% as expected

clearvars;
close all;
%part one run with nominal values
w = 6.5*(2 * pi /60);
r2 = 0.8;
a1 = 0.058;
T = 10000;
r1 = 4.3;
fun = @(t,y) ode_(w, r2, a1, y, t, r1);
Tau = T;
[t,y_2] = ode45(fun, [0,Tau],[0,0]);
y = y_2(:,1);
y_m = 1/Tau * trapz(t,y);
st = 3* w * (1/Tau * trapz(t,(y - y_m).^2)).^(1/2);
% Plot the ODE to ensure it looks chaotic
xlim([100, 200])
plot(t,y)

% part 2: run a sweep of omega
r2 = 0.8;
a1 = 0.058;
T = 10000;
r1 = 4.3;
% change w_lst to change what range to sweep 
w_lst = 1:10:100; % .3:.001:.9
st_ = zeros(length(w_lst),1);
i =1;
for w = w_lst
    fun = @(t,y) ode_(w, r2, a1, y, t, r1);
    Tau = 3 * w * T;
    tspin = .1 * Tau;
    [t,y_2] = ode45(fun, [tspin,Tau],[0,0]);
    y = y_2(:,1);
    y_m = 1/Tau * trapz(t,y);
    st_(i) = 3* w * (1/Tau * trapz(t,(y - y_m).^2)).^(1/2);
    i = i+1;
end
plot(w_lst, st_)
xlabel('Omega (rad/sec)')
ylabel('Standard Devation(rad/sec)')
