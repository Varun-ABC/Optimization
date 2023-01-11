clearvars 
close all

Nelem = 15;
L = 7.5; % semi-span in meters
x = [0:L/Nelem:L].';
r_in = ones(Nelem+1,1) * .0415;
r_out = ones(Nelem+1,1) * .05;
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