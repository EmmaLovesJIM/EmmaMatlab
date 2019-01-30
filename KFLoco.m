%% Prepare
clear all
close all
clc
m = 1800;
N = 1000;
start = 25;
sigma = 0.5;

%% System model
sys = tf([1/m], [1 0.001 0]);
sysd = c2d(sys, 0.1);
sysSS = ss(sysd);
sysSS.C(2) = 0;

%% System simulation
u=[-1*ones(start,1);ones(N-start,1)];
u = 10*sin(.01*(1:N));

noise=randn(N,1); noise=sigma*(noise-mean(noise));

[y, t, x] = lsim(sysSS,u);

y0 = y;
y = y + noise;

%% Set up for Kalman Filter
A = sysSS.A;
B = sysSS.B;
C = sysSS.C';
n = length(sysSS.A);
% Set up the matrices and vectors for Kalman filter
x = zeros(n,1);
P = eye(n)*1e-3;
% For parameter estimation
n2 = length(sysSS.A) + length(sysSS.B);
P2 = eye(n2)*1e-3;
Rw2 = eye(n2);
x2 = zeros(n2,1);
theta = x2;
iae =0;
% Noise information
rv = 1;  Rw = 0.01*eye(n); xSave=[];

for t= 3:N
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
    % KF for parameter estimation
    yt = sysSS.C*x;

%	Kalman filter algorithm
%	Prediction
	P2 = P2 + Rw2;
% 	Correction
	K2 = P2*x2/(rv +x2'*P2*x2);
	perror = (yt - x2'*theta);
	theta = theta + K2*perror;
	P2 = (eye(n2) - K2*x2')*P2;

%	Store the predicted error
	iae = iae + abs(perror);

% 	Time shift the variables
	yt_1 = yt;
	ut_1 = u(t);

% 	Update the observation vector
	x2(2) = x2(1); x2(1) = -yt_1;
	x2(4) = x2(3); x2(3) = ut_1; 

%	Store the output of the system
	store(t,1) = yt;

%	Store the system parameters
	store(t,2:n2+1) = theta';
end

L = 2;

yest = sysSS.C*xSave;

h = figure(1);
plot([1:N],y, 'LineWidth', L)
hold on
plot([3:N], yest,'LineWidth', L)
plot([1:N],y0, 'LineWidth', L)
%plot([1:N],u, 'LineWidth', L)
xlabel('Iterations')
ylabel('System response')
legend('y_{meas}', 'y_{est}', 'y_{true}','u');
FS = findall(h, '-property', 'Fontsize');
set(FS, 'FontSize', 14);