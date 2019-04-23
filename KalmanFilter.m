function adapt_lms(block)
% Level-2 MATLAB file S-Function for system identification using 
% Least Mean Squares (LMS).

%   Copyright 1990-2011 The MathWorks, Inc.

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register dialog parameter: LMS step size 
  block.NumDialogPrms = 7;
  block.DialogPrmsTunable = {'Tunable', 'Tunable', 'Tunable', 'Tunable', 'Tunable', 'Tunable', 'Tunable'};
  %block.DialogPrm(1).Name = 'A';
  block.DialogPrm(1).DataTypeId = 0;
  %block.DialogPrm(2).Name = 'B';
  block.DialogPrm(2).DataTypeId = 0;
  %block.DialogPrm(3).Name = 'C';
  block.DialogPrm(3).DataTypeId = 0;
  %block.DialogPrm(4).Name = 'P';
  block.DialogPrm(4).DataTypeId = 0;
  %block.DialogPrm(5).Name = 'rv';
  block.DialogPrm(5).DataTypeId = 0;
  %block.DialogPrm(6).Name = 'Rw';
  block.DialogPrm(6).DataTypeId = 0;
  block.DialogPrm(7).DataTypeId = 0;
  
  %% Regieste number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;

  block.InputPort(1).Complexity   = 'Real'; 
  block.InputPort(1).DataTypeId   = 0;
  block.InputPort(1).SamplingMode = 'Sample';
  block.InputPort(1).Dimensions   = 1;
  
  block.InputPort(2).Complexity   = 'Real';
  block.InputPort(2).DataTypeId   = 0;
  block.InputPort(2).SamplingMode = 'Sample';
  block.InputPort(2).Dimensions   = 1;
  
  block.OutputPort(1).Complexity   = 'Real';
  block.OutputPort(1).DataTypeId   = 0;
  block.OutputPort(1).SamplingMode = 'Sample';
  block.OutputPort(1).Dimensions   = 2;
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  % block.RegBlockMethod('CheckParameters',         @CheckPrms);
  block.RegBlockMethod('ProcessParameters',       @ProcessPrms);
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('Start',                   @Start);  
  block.RegBlockMethod('WriteRTW',                @WriteRTW);
  block.RegBlockMethod('Outputs',                 @Outputs);
  block.RegBlockMethod('Terminate', @Terminate); % Required

  
%endfunction

function DoPostPropSetup(block)

  %% Setup Dwork  
  block.NumDworks = 2;
  n = 2; % System dimension
  block.Dwork(1).Name            = 'x';
  block.Dwork(1).Dimensions      = n;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;
  
  block.Dwork(2).Name            = 'P';
  block.Dwork(2).Dimensions      = n*n;
  block.Dwork(2).DatatypeID      = 0;      % double
  block.Dwork(2).Complexity      = 'Real'; % real
  block.Dwork(2).UsedAsDiscState = true;

  %% Register all tunable parameters as runtime parameters.
  block.AutoRegRuntimePrms;

%endfunction

function Start(block)
  
  %% Initialize Dwork 
  block.Dwork(1).Data = [0, block.DialogPrm(7).Data];
  block.Dwork(2).Data = [block.DialogPrm(4).Data(:,1); block.DialogPrm(4).Data(:,2)];
  
%endfunction

function Update(block)

%endfunction

function Outputs(block)
A = block.DialogPrm(1).Data;
B = block.DialogPrm(2).Data;
C = block.DialogPrm(3).Data;
rv = block.DialogPrm(5).Data;
Rw = block.DialogPrm(6).Data;

x = block.Dwork(1).Data;
P = [block.Dwork(2).Data(1:2), block.Dwork(2).Data(3:4)];
u = block.InputPort(1).Data;
%disp(['u: ', num2str(u)])
y = block.InputPort(2).Data;
%disp(['y: ', num2str(y)])
%% Kalman filter for state estimation
%	Prediction
	x = A*x + B*u; % x^(t+1|t)= P x(t|t) + Q u(t)
    % Rw process noise covariance matrix
 	P = A*P*A' + Rw; % phi(t+1|t) = P phi(t|t) P' + Rw
%	Correction
    % rv ouput noise variance
	K = P*C'/(rv + C*P*C'); % K(t+1) = [phi(t+1|t) H' ] / [rv + H phi(t+1|t) H']
	x = x + K*(y - C*x); % x^(t+1|t+1) = x^(t+1|t) + K(t+1) [y(t+1) - H x^(t+1|t)]
	P = (eye(2) - K*C)*P; % phi(t+1|t+1) = [I - K(t) H] phi(t+1|t)

block.Dwork(2).Data = [P(:,1); P(:,2)];
block.Dwork(1).Data = x;
block.OutputPort(1).Data = x;
  
  
%endfunction

function Terminate(block)

%end Terminate
