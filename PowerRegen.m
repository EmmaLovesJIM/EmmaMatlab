clear all, close all, clc
mL = 650;
mW = 1200;
 
dt = 0.1;
a = 0.3;
t = 0;
i = 1;
v = 15/3.6;

while v > 0
    v = v-a*dt;
    P(i) = (1.1*mL+1.04*mW)*v*a;
    T(i) = t;
    t = t+dt;
    i = i+1;
end

plot(T, P)