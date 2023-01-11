% this is test code to ensure that the standard devation matches the
% nominal
close all 
clearvars
max_p = 6;
times = zeros(1, max_p);
stresses = zeros(1, max_p);
sds = zeros(1, max_p);

for n = 1:max_p
    r_in = ones(Nelem+1,1) * .0415;
    r_out = ones(Nelem+1,1) * .05;
    [stress,sd] = standard_dev(n, r_in, r_out);
    sds(n) = sd(1);
    f = @() standard_dev(n, r_in, r_out);
    times(n) = timeit(f);
end
percent_error = [0,0,0,0,0];
for i = 2:max_p
percent_error(i-1) = abs(sds(6) - sds(i)) / sds(6) * 100;
end
figure 
plot([1:max_p], times)
title('Effect of number of points on computational time')
xlabel('Number of quad. points')
ylabel('Time in seconds')
figure
plot([2:max_p], percent_error)
title('Effect of number of points on accuracy of standard devation')
xlabel('Number of quad. points')
ylabel('Percent error')
