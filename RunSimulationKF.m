%% Initialise
clear all, close all, clc
mL = 650;
mWz = 1200;
v0 = 10/3.6;
s0 = 0;
Ts = 0.01;
Pres = 1/1000*[5.7/771 0 1.6]; %Strahl formula for m/s velocity
i = -5/1000; %Gradient, uphill positive
%System model
%sys = tf([1/(mL+mWz)], [1 0.01/(mL+mWz) 0]);
%sysSSc = ss(sys);
A = [0 1; 0 -0.06/(mL+mWz)]; B = [0;1/(mL+mWz)]; C = [1 0];
sysSSc = ss(A, B, C, []);
sysSS = c2d(sysSSc, Ts);
sysSS.C = [0,1];
%Ts = sysSS.Ts;
n = 2;
P = eye(n)*1e3;
rv = 20;
smax = 20;
Kp = 600;
Ki = 20;
P1 = 1e-1;
P2 = 0;
%Rw = 1*[P2 0; 0 P1];
Rw = [1e-3 0; 0, 1e-3];
P = eye(n)*1e-3;

%% Run Simulation
tmax = 20;
nmax = 200;
t = linspace(0, tmax, nmax);
u = 300*idinput(nmax);
simin.time = t;
simin.signals.values = [-300*ones(nmax,1)];

sim('SimulationKF.slx')
%% Plot
L = 2;

xest = stateout.Data(1:end,1);
vest = stateout.Data(1:end,2);
N = length(xest);
tplot = linspace(0, tmax, N);
y = simout.Data(:,2);
%mest = thetaSave.^(-1);

figure
subplot(5,1,1)
plot(simout.Time,simout.Data(:,2), 'LineWidth', L);
hold on
plot(tplot, vest,'LineWidth', L)
plot(simout.Time,simout.Data(:,4), 'LineWidth', L);
ylim([0 ceil(v0)])
%plot(tplot, vest, 'LineWidth', L)
%xlabel('t/s')
ylabel('Velocity')
legend('v_{true}', 'v_{est}', 'v_{obs}');
grid on
subplot(5,1,2)
%plot(tplot,y, 'LineWidth', L);
hold on
plot(tplot, xest,'LineWidth', L)
%plot(tplot,x0, 'LineWidth', L)
plot(simout.Time, simout.Data(:,3), 'LineWidth', L)
xlabel('t/s')
ylabel('Position')
legend('x_{est}', 'x_{true}');
grid on
subplot(5,1,3:4)
%plot(tplot,y, 'LineWidth', L);
hold on
plot(xest, vest,'LineWidth', L)
%plot(tplot,x0, 'LineWidth', L)
plot(simout.Data(:,3), simout.Data(:,2), 'LineWidth', L)
xlabel('s/m')
ylabel('v/m/s')
ylim([0 ceil(v0)])
legend('bc_{est}', 'bc_{true}');
title(['Braking distance ', num2str(max(simout.Data(:,3))), ' m'])
grid on
subplot(5,1,5)
plot(simout.Data(:,3), simout.Data(:,5)/(-10*mL), 'LineWidth', L)
ylabel('\mu')
xlabel('s/m')



