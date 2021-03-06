%% Initialise
clear all, close all, clc
mL = 650;
mWz = 1200;
v0 = 12/3.6;
s0 = 0;
Ts = 0.010;
Pres = 1/1000*[5.7/771 0 1.6]; %Strahl formula for m/s velocity
i = -5/1000; %Gradient, uphill positive
ssigma = 1e-1; % Slip at 68% of mumax, e.g. dry: 1%, wet: 10%
%System model
%sys = tf([1/(mL+mWz)], [1 0.01/(mL+mWz) 0]);
%sysSSc = ss(sys);
A = [0 1; 0 -0.06/(mL+mWz)]; B = [0;1/(mL+mWz)]; C = [1 0];
sysSSc = ss(A, B, C, []);
sysSS = c2d(sysSSc, Ts);
sysSS.C = [0,1];
%Ts = sysSS.Ts
n = 2;


% Braking parameters
mumax = 0.15;
smax = 25;
soffset = 0;
% Controller parameters
%pmax: 3500
Kp = 0.2*3500;
Ki = 50;
% Noise parameters
P1 =0e-2;
P2 = 0;
% Kalman filter parameters
Rw = [1e-3 0; 0, 1e-3];
P = eye(n)*1e-3;
rv = 20;
% Pre-braking parameters
alpha = 2;
beta = -.2;

%% Run Simulation
tmax = 25;
nmax = 200;
t = linspace(0, tmax, nmax);
u = 300*idinput(nmax);
simin.time = t;
simin.signals.values = [-300*ones(nmax,1)];

sim('SimulationKFslip.slx')
%% Plot
L = 2;

xest = stateout.Data(1:end,1);
vest = stateout.Data(1:end,2);
N = length(xest);
tplot = linspace(0, tmax, N);
y = simout.Data(:,2);
%mest = thetaSave.^(-1);

figure
subplot(3,2,1)
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
subplot(3,2,2)
%plot(tplot,y, 'LineWidth', L);
hold on
plot(tplot, xest,'LineWidth', L)
%plot(tplot,x0, 'LineWidth', L)
plot(simout.Time, simout.Data(:,3), 'LineWidth', L)
xlabel('t/s')
ylabel('Position')
legend('x_{est}', 'x_{true}');
grid on
subplot(3,2,3)
plot(simout.Data(:,3), simout.Data(:,5)/(-10*mL), 'LineWidth', L)
ylabel('\mu')
xlabel('s/m')
ylim([0 0.25])
grid on
subplot(3,2,[4,6])
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
subplot(3,2,5)
plot(simout.Data(:,3), simout.Data(:,[1,6:7]), 'LineWidth', L)
ylabel('a')
xlabel('s/m')
legend('a_{true}', 'a_{est}', 'a_{req}');
ylim([-1,0])
grid on



