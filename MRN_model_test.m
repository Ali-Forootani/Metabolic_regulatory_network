clear all; close all; clc; 
figpath = '../figures/';
addpath('./utils');
addpath('./bioutils');
changeplot;

plottag = 2;

%% generate Data

% Size of system
n = 8;  % Update to match the number of state variables in MRN model

% Integrate
dt = 0.1;
tspan = 0:dt:5; % time vector
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% Define initial conditions
rng(1); % seed the random number generator


%%
% sss = 300; % 3*sss is the number of initial conditions
%Sinit = rand(sss,n);
%Sinit = [Sinit; 2*rand(sss,n)];
%Sinit = [Sinit; 3*rand(sss,n)];
% Sinit= lhsdesign(sss,n);  %%%% latin hybercube sample
%Sinit;
%%
% Define initial conditions to match Python code
C1_0 = 500;   % Initial concentration for C1
C2_0 = 1000;  % Initial concentration for C2
M_0 = 20;      % Initial concentration for M
Q_0 = 0.1;    % Initial concentration for Q
R_0 = 0.01;    % Initial concentration for R
T1_0 = 0.05;  % Initial concentration for T1
T2_0 = 0.05;  % Initial concentration for T2
RP_0 = 0.05;   % Initial concentration for RP

Sinit = [C1_0, C2_0, M_0, Q_0, R_0, T1_0, T2_0, RP_0];
 sss = 300; % 3*sss is the number of initial conditions
% Sinit = rand(sss,n);
Sinit = [Sinit; 2*rand(sss,n)];
Sinit = [Sinit; 3*rand(sss,n)];
%Sinit= lhsdesign(sss,n);  %%%% latin hybercube sample
Sinit;
%%
% 
% % Define base initial conditions to match Python code
% C1_0 = 500;   % Initial concentration for C1
% C2_0 = 1000;   % Initial concentration for C2
% M_0 = 1;    % Initial concentration for M
% Q_0 = 0.2;  % Initial concentration for Q
% R_0 = 0.1;  % Initial concentration for R
% T1_0 = 0.01; % Initial concentration for T1
% T2_0 = 0.05; % Initial concentration for T2
% RP_0 = 0.1;  % Initial concentration for RP
% 
% % Base initial conditions
% y0 = [C1_0, C2_0, M_0, Q_0, R_0, T1_0, T2_0, RP_0];
% 
% % Number of different initial conditions to generate
% num_conditions = 300;
% 
% % Define perturbation range (small noise)
% perturbation = 0.02; % 10% variation around the initial value
% 
% % Generate initial conditions with small perturbations
% Sinit = zeros(num_conditions, length(y0)); % Preallocate
% 
% for i = 1:num_conditions
%     Sinit(i, :) = y0 + perturbation * y0 .* (2*rand(1, length(y0)) - 1)
% end
% 
% Sinit;  % Display generated initial conditions

%%

% % Define base initial conditions to match Python code
% C1_0 = 500;   % Initial concentration for C1
% C2_0 = 1000;  % Initial concentration for C2
% M_0 = 20;      % Initial concentration for M
% Q_0 = 0.1;    % Initial concentration for Q
% R_0 = 0.01;    % Initial concentration for R
% T1_0 = 0.05;  % Initial concentration for T1
% T2_0 = 0.05;  % Initial concentration for T2
% RP_0 = 0.05;   % Initial concentration for RP
% 
% 
% % Base initial conditions
% y0 = [C1_0, C2_0, M_0, Q_0, R_0, T1_0, T2_0, RP_0];
% 
% % Number of different initial conditions to generate
% num_conditions = 200;
% 
% % Define individual perturbation ranges for each variable
% perturbations = [2, 5, 0.05, 0.1, 0.1, 0.05, 0.05, 0.01];  % Each value corresponds to a variable in y0
% 
% % Generate initial conditions with individual perturbations
% Sinit = zeros(num_conditions, length(y0)); % Preallocate
% 
% for i = 1:num_conditions
%     % Apply different perturbation for each variable
%     Sinit(i, :) = y0 + perturbations .* y0 .* (2*rand(1, length(y0)) - 1);
% end
% 
% Sinit;  % Display generated initial conditions



%%
measure = size(Sinit);

for ii = 1:measure(1)
    % Integrate for each initial condition using MRN model
    [t1,x1] = ode45(@MRN_Model, tspan, Sinit(ii,:));
    % Store each instance
    tt(:,ii) = t1;
    x(:,:,ii) = x1;
end

%% Add noise

eps = 0; % no error used here
xn = x + eps*randn(size(x));

% Compute Derivative
% Calculate exactly and add error
xt = []; dxt= []; t = [];

for ll = 1:measure(1)
    for ii = 1:length(tspan)
        dxf(ii,:,ll) = MRN_Model(tspan(ii), xn(ii,:,ll)); % Use tspan(ii) for correct time input
    end
    epsdt = 0;
    dxf = dxf + epsdt*randn(size(dxf));
    
    dxt = [dxt; dxf(:,:,ll)];
    xt = [xt; xn(:,:,ll)];
    t = [t; tt(:, ll)];
end

%% Adjust num2plot dynamically
num2plot = 1:min(5*length(tspan), length(t));

if plottag > 1
    % Plot data and derivatives
    figure(2)
    plot(t(num2plot), xt(num2plot,:),'o')  % Changed t(num2plot,:) to t(num2plot)
    hold on
    xlabel('time')
    ylabel('concentrations')
    drawnow
    
    figure(6)
    plot(t(num2plot), dxt(num2plot,:), '.')  % Changed t(num2plot,:) to t(num2plot)
    hold on
    xlabel('time')
    ylabel('derivative of concentrations w/ time')
    drawnow
end
