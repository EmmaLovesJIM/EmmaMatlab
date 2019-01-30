%% Initialise
clear all, close all, clc
mL = 650;
mWz = 1200;
v0 = 10/3.6;
s0 = 0;
Pres = 1/1000*[5.7/771 0 1.6]; %Strahl formula for m/s velocity
i = -5/1000; %Gradient, uphill positive
%System model
sys = tf([1/(mL+mWz)], [1 0.01/(mL+mWz) 0]);
sysd = c2d(sys, 0.1);
sysSS = ss(sysd);
Ts = sysd.Ts;



%% Run Simulation
tmax = 21;
nmax = 200;
t = linspace(0, tmax, nmax);
simin.time = t;
simin.signals.values = -300*ones(nmax,1);

sim('Simulation.slx')
%% Set up for Kalman Filter
u = simout.Data(:,5);
y0 = simout.Data(:,3);
% y = y0+randn(size(y0));
% A = sysSS.A;
% B = sysSS.B;
% C = sysSS.C';
% n = length(sysSS.A);
% N = length(y);
y = y0+randn(size(y0));
A = [0, -1;1, 2];
B = 2.703e-6*ones(2,1);
C = [0 1]';
n = length(sysSS.A);
N = length(y);
% Set up the matrices and vectors for Kalman filter
x = [0;0];
P = eye(n)*1e3;
% Noise information
rv = 10;  Rw = 1e-6*eye(n); xSave=[];

for t= 2:N
%	Kalman filter for state estimation
%	Prediction
	x = A*x + B*u(t-1); % x^(t+1|t)= P x(t|t) + Q u(t)
    % Rw process noise covariance matrix
 	P = A*P*A' + Rw; % phi(t+1|t) = P phi(t|t) P' + Rw
%	Correction
    % rv ouput noise variance
	K = P*C/(rv + C'*P*C); % K(t+1) = [phi(t+1|t) H' ] / [rv + H phi(t+1|t) H']
	x = x + K*(y(t) - C'*x); % x^(t+1|t+1) = x^(t+1|t) + K(t+1) [y(t+1) - H x^(t+1|t)]
	P = (eye(n) - K*C')*P; % phi(t+1|t+1) = [I - K(t) H] phi(t+1|t)
    xSave=[xSave, x];
end

L = 2;

yest = C'*xSave;
tplot = linspace(0, tmax, N);

figure
plot(tplot,y, 'LineWidth', L);
hold on
plot(tplot(2:end), yest,'LineWidth', L)
plot(tplot,y0, 'LineWidth', L)
%plot([2:N],sysSS.C.*xSave', 'LineWidth', L)
%plot(simout.time, sysSS.C*simout.Data(:, 7:8)', 'LineWidth', L)
plot(tplot(3:N), diff(yest), 'LineWidth', L)
xlabel('Iterations')
ylabel('System response')
legend('y_{meas}', 'y_{est}', 'y_{true}','Simulink');


%% Plot
figure
subplot(321)
plot(simout.time, simout.Data(:,1))
ylabel('a')
%
subplot(323)
plot(simout.time, simout.Data(:,2))
ylabel('v')
%
subplot(325)
plot(simout.time, simout.Data(:,3))
ylabel('s')
xlabel('t')
%
subplot(322)
plot(simout.time, simout.Data(:,5:6))
legend('F_res', 'F_t/db')
ylabel('F')
%
subplot(324)
plot(simout.time, simout.Data(:,4))
ylabel('F_{ww}')
xlabel('t')
%
subplot(326)
plot(simout.Data(:,3), simout.Data(:,2))
ylabel('v')
xlabel('s')