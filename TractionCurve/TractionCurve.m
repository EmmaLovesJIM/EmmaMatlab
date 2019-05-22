clear all, close all, %clc

I = 57.9; %A
U = 50; %V

mu = 0.4;
m = 600;
mW = 1200;
Fmax = mu*m*10;

P = U*I;

v = linspace(0, 15/3.6);

F = P./v;
Fres = 10*(m+mW)*1/1000*(5.7/771*(3.6*v).^2 +1.6)+10*(m+mW)*[0.01;0.02];

Ind = find(F > Fmax);
F(end) = 0;

F(Ind) = Fmax;

L = 2;
plot(3.6*v, F, 'LineWidth', L)
hold on;
plot(3.6*v, Fres, 'LineWidth', L)
ts = ['Traction curve, $P =$' num2str(U*I) ' W, $\mu =$' num2str(mu)  ', $m_L =$' num2str(m) ' kg, $m_{Wz} =$' num2str(mW) ' kg'];
title(ts ,'interpreter','latex')
legend('F_{t/b}', 'F_{R, 1%}', 'F_{R, 2%}')
xlim([0 16])
xlabel('$v /\left(kmh^{-1}\right)$','interpreter','latex')
ylim([0 1.1*Fmax])
ylabel('$F$/N','interpreter','latex')
grid on
