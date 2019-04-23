%% Initialise
clear all, close all, clc
mL = 650;
mWz = 1200;
v0 = 10/3.6;
s0 = 0;
Ts = 0.01;
Pres = 1/1000*[5.7/771 0 1.6]; %Strahl formula for m/s velocity
i = 0%-5/1000; %Gradient, uphill positive
%System model
%sys = tf([1/(mL+mWz)], [1 0.01/(mL+mWz) 0]);
%sysSSc = ss(sys);
A = [0 1; 0 -0.06/(mL+mWz)]; B = [0;1/(mL+mWz)]; C = [1 0];
sysSSc = ss(A, B, C, []);
sysSS = c2d(sysSSc, Ts);
%Ts = sysSS.Ts;



%% Run Simulation
tmax = 20;
nmax = 200;
t = linspace(0, tmax, nmax);
u = 300*idinput(nmax);
simin.time = t;
simin.signals.values = [-300*ones(nmax,1)];

sim('Simulation.slx')
%% Set up for Kalman Filter
sigma = 0; % Noise variance
u = simoutD.Data(:,5);
y0 = simoutD.Data(:,2);
y = y0+sigma*randn(size(y0)); %Systematic offset +.5*linspace(v0,0,length(u))';
x0 = simoutD.Data(:,3);
A = sysSS.A;
B = sysSS.B;
C = sysSS.C';
C = [0;1];
n = length(sysSS.A);
N = length(y);
% Set up the matrices and vectors for Kalman filter
x = [0;v0];
P = eye(n)*1e-3;
% Noise information
rv = 0.1;  Rw = 1e-2*eye(n); xSave=[0;0];

% RLS setup
nRLS = 1;
xRLS = zeros(nRLS,1);
theta = [1/(mL+mWz)];
thetaSave = [];

% Set up the covariance matrix
PRLS = 1e8*eye(nRLS);

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
    
    % RLS for parameter estimation 
    if abs(u) > 0
        a = xSave(1,t-1) - xSave(2,t);
    % Observation vector
        xRLS(1) = u(t-1)/a;
    
    
        PRLS = PRLS - PRLS*xRLS*((1 + xRLS'*PRLS*xRLS)^(-1))*xRLS'*PRLS;
        theta = theta + PRLS*xRLS*(y(t) - xRLS'*theta);
    end
    thetaSave = [thetaSave, theta];
 
end

L = 2;

xest = xSave(1,2:end);
vest = xSave(2,2:end);
tplot = linspace(0, tmax, N);
mest = thetaSave.^(-1);

figure
subplot(5,1,1:2)
plot(tplot,y, 'LineWidth', L);
hold on
plot(tplot(2:end), vest,'LineWidth', L)
plot(tplot,y0, 'LineWidth', L)
plot(tplot(2:end), vest - y0(2:end)', 'LineWidth', L)
%xlabel('t/s')
ylabel('Velocity')
legend('v_{meas}', 'v_{est}', 'v_{true}', 'Delta');
grid on
subplot(5,1,3:4)
%plot(tplot,y, 'LineWidth', L);
hold on
plot(tplot(2:end), xest,'LineWidth', L)
plot(tplot,x0, 'LineWidth', L)
plot(tplot(2:end), xest - x0(2:end)', 'LineWidth', L)
xlabel('t/s')
ylabel('Position')
legend('x_{est}', 'x_{true}', 'Delta');
grid on
subplot(515)
plot(tplot(2:end), -1*mest, 'LineWidth', L)
hold on
plot(tplot([1,end]), [mL+mWz, mL+mWz], 'LineWidth', L)
ylim([0,5000])
legend('m_{est}', 'm_{true}')
grid on



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
legend('F_{res}', 'F_{t/db}')
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