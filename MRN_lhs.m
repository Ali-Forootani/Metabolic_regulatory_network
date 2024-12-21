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
tspan = 0:dt:15; % time vector
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
T1_0 = 0.01;  % Initial concentration for T1
T2_0 = 0.01;  % Initial concentration for T2
RP_0 = 0.01;   % Initial concentration for RP

%Sinit = [C1_0, C2_0, M_0, Q_0, R_0, T1_0, T2_0, RP_0];
% sss = 300; % 3*sss is the number of initial conditions
%  Sinit = rand(sss,n);
% Sinit = [Sinit; 2*rand(sss,n)];
%  Sinit = [Sinit; 3*rand(sss,n)];
% %Sinit= lhsdesign(sss,n);  %%%% latin hybercube sample
%  Sinit;

%%%%%% latin hybercube sample

nSamples = 1000;    % Total number of design points (samples)
nVariables = 8;   % Number of variables

% Initial condition (specific starting point)
Sinit = [C1_0, C2_0, M_0, Q_0, R_0, T1_0, T2_0, RP_0];  % Your predefined initial condition for each variable

% Generate Latin Hypercube Sample (LHS) for the remaining points
lhs_rest = lhsdesign(nSamples - 1, nVariables);

% Specify the bounds for scaling  ORIGINAL one this good for s1,s2,s4,s5,s6
% lowerBounds = [400, 900, 19, 0.05, 0   , 0,  0.01, 0.01];  % Lower bound for each variable for s4 it can be 0.05 and up 0.2
% upperBounds = [600, 1100, 25, 0.2, 0.1, 0.1, 0.1, 0.1];     % Upper bound for each variable

%%%%%
lowerBounds = [400, 900, 18, 0.05, 0  ,  0,  0.05, 0.01];  % Lower bound for each variable for s4 it can be 0.05 and up 0.2
upperBounds = [600, 1100, 22, 0.2, 0.1,  0.1,  0.8,    0.1];     % Upper bound for each variable
% Scale the LHS points for the rest of the design
lhs_rest_scaled = bsxfun(@plus, lowerBounds, bsxfun(@times, lhs_rest, (upperBounds - lowerBounds)));

% Combine initial condition as the first row with the rest of the scaled LHS points
Sinit= [Sinit; lhs_rest_scaled];





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
    % Define font sizes and line width
    fontSizeLabel = 12;
    fontSizeTitle = 14;
    fontWeight = 'bold';
    lineWidth = 2;

    % Plot data and derivatives
    figure(2)
    plot(t(num2plot), xt(num2plot, :), 'o', 'LineWidth', lineWidth)  % Changed t(num2plot,:) to t(num2plot)
    hold on
    xlabel('Time (min)', 'FontSize', fontSizeLabel, 'FontWeight', fontWeight)
    ylabel('Concentrations', 'FontSize', fontSizeLabel, 'FontWeight', fontWeight)
    title('Concentration over Time', 'FontSize', fontSizeTitle, 'FontWeight', fontWeight)
    set(gca, 'FontSize', 10)  % Set axis tick font size
    grid on
    drawnow

    figure(6)
    plot(t(num2plot), dxt(num2plot, :), '.', 'LineWidth', lineWidth)  % Changed t(num2plot,:) to t(num2plot)
    hold on
    xlabel('Time (min)', 'FontSize', fontSizeLabel, 'FontWeight', fontWeight)
    ylabel('Derivative of Concentrations w/ Time', 'FontSize', fontSizeLabel, 'FontWeight', fontWeight)
    title('Derivative of Concentrations over Time', 'FontSize', fontSizeTitle, 'FontWeight', fontWeight)
    set(gca, 'FontSize', 10)  % Set axis tick font size
    grid on
    drawnow
end


