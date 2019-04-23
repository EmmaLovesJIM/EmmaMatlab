% EKF estimates the state of a system using the 
%     extended Kalman filter
%

% Clear the previous set of parameters
clear;

% System parameters
a1 = -1.5;
a2 = 0.7;
b0 = 0.9;

% Zero previous input and output variables
yt_1 = 0; yt_2 = 0;
ut_1 = 0; 

% Input to the system is a step
ut = 1.0; start = 25;

% Set up the noise
randn('seed',0);


% Set up cbar
cbar = [0 1 0 0 0]';

% Set up initial vectors
n = 3 + 2;
z = zeros(n,1);

% Set up the covariance matrix
P = 10000*eye(n);

% Set up Rw matrix
rv = 1;
Rw = eye(n);

for i= 1:50

%	Change the reference signal
	if i == start
	  ut = -1;
	end

%	System model
	noise = 0.25*rand;
	yt = -a1*yt_1 -a2*yt_2 + b0*ut_1 + noise;

%	Update Abar, bbar and Jacobian
	Abar = [0 -z(4) 0 0 0
		1 -z(3) 0 0 0
		0 0     1 0 0
 		0 0     0 1 0
        	0 0     0 0 1];
	bbar = [0 z(5) 0 0 0]';
	J = [0 -z(4) 0      -z(2) 0
	     1 -z(3) -z(2)  0     ut_1
	     0 0     1      0     0
 	     0 0     0      1     0
             0 0     0      0     1];

%	Extended Kalman filter
%	Prediction
	z = Abar*z + bbar*ut_1;
	P = J*P*J' + Rw;
% 	Correction
	K = P*cbar/(rv + cbar'*P*cbar);
	perror = (yt - cbar'*z);
	z = z + K*perror;
	P = (eye(n) - K*cbar')*P;

% 	Time shift the variables
	yt_2 = yt_1;
	yt_1 = yt;
	ut_1 = ut;

%	Store the output of the system
	store(i,1) = yt;
	store(i,2) = ut;
	store(i,3) = z(2);
end

% Clear the graphics screen
clf;
t = 0:1:i-1;
% Determine the noise free estimate
yf = filter([0 b0],[1 a1 a2],store(:,2));
plot(t,store(:,1),t,store(:,3),'o',t,yf,'*');
xlabel('Iterations'),ylabel('System response');

