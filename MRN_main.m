clear all; close all; clc; 
figpath = '../figures/';
addpath('./utils');
addpath('./bioutils');
changeplot;

plottag = 2;

%% generate Data

%size of system
n = 8;  % Update to match the number of state variables in MRN model

% Integrate
dt = 0.1;
tspan = 0:dt:5; % time vector
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% Define initial conditions
rng(1); % seed the random number generator

sss = 300; % 3*sss is the number of initial conditions
Sinit = rand(sss,n);
Sinit = [Sinit; 2*rand(sss,n)];
Sinit = [Sinit; 3*rand(sss,n)];
%Sinit= lhsdesign(sss,n);  %%%% latin hybercube sample

measure = length(Sinit);

for ii = 1:measure
    % Integrate for each initial condition using MRN model
    [t1,x1] = ode45(@MRN_Model, tspan, Sinit(ii,:));
    % Store each instance
    tt(:,ii) = t1;
    x(:,:,ii) = x1;
end

%% add noise

eps = 0; % no error used here
xn = x + eps*randn(size(x));

% Compute Derivative
% Calculate exactly and add error
xt = []; dxt= []; t = [];

for ll = 1:measure
    for ii = 1:length(tspan)
        dxf(ii,:,ll) = MRN_Model(t, xn(ii,:,ll));
    end
    epsdt = 0;
    dxf = dxf + epsdt*randn(size(dxf));
    
    dxt = [dxt; dxf(:,:,ll)];
    xt = [xt; xn(:,:,ll)];
    t = [t; tt(:, ll)];
end

num2plot = 1:5*length(tspan);

if plottag > 1
    % Plot data and derivatives
    figure(2)
    plot(t(num2plot,:) ,xt(num2plot,:),'o')
    hold on
    xlabel('time')
    ylabel('concentrations')
    drawnow
    
    figure(6)
    plot(t(num2plot,:), dxt(num2plot,:), '.')
    hold on
    xlabel('time')
    ylabel('derivative of concentrations w/ time')
    drawnow
end
