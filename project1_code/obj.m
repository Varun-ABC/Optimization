function [r] = obj(a, h, l_iphone,w_ip, Lg1, Lg2, kg, kal)
rtot = resistance(h, l_iphone,w_ip, Lg1, Lg2, kg, kal, a(1), a(2), a(3), a(4));
r  = 1/rtot;