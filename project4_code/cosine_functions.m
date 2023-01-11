clearvars
close all

L = 7.5; % semi-span in meters
x = [0:L/100:L].';
c1 = cos((2-1)* pi * x / (2*L));
c2 = cos((2*2-1)* pi * x / (2*L));
c3 = cos((2*3-1)* pi * x / (2*L));
c4 = cos((2*4-1)* pi * x / (2*L));
c = [c1,c2,c3,c4];
figure
plot(x, c(:,1))
title('Cosine function when n=1')
xlabel('Distance along span (m)') 
ylabel('Unscaled peturabtion') 
figure
plot(x, c(:,2))
title('Cosine function when n=2')
xlabel('Distance along span (m)') 
ylabel('Unscaled peturabtion') 
figure
plot(x, c(:,3))
title('Cosine function when n=3')
xlabel('Distance along span (m)') 
ylabel('Unscaled peturabtion') 
figure
plot(x, c(:,4))
title('Cosine function when n=4')
xlabel('Distance along span (m)') 
ylabel('Unscaled peturabtion') 