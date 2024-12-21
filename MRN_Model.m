function dS = MRN_Model(t, S)

% Unpack state variables
C1 = S(1); % Species 1
C2 = S(2); % Species 2
M = S(3);  % Intermediate Species
Q = S(4);  % Species Q
R = S(5);  % Species R
T1 = S(6); % Species T1
T2 = S(7); % Species T2
RP = S(8); % Regulatory Protein RP

% Updated Constants
kcat1 = 3000;  % mM/min
kcat2 = 2000;  % mM/min
kr = 1260;        % mM/min
KT = 1000;     % mM
Kr = 7;     % mM
kde = 0.01;    % 1/min
kdRP = 0.2;    % 1/min
kdT1 = 0.05;   % 1/min
kdT2 = 0.05;   % 1/min
nQ = 300;      % mM
nRP = 300;     % mM
nT1 = 400;     % mM
nT2 = 1500;    % mM
nR = 7459;     % mM
beta_p = 1;    % assumed
gamma = 20;    % mM
alpha = 0.03;  % mM
%%

%%
 %%%%  Hill function
%  h=2;
%  vRP = C1.^h./(20^h+C1.^h);        %%%% \gamma
%  vT_2 = 0.03^h./(0.03^h+RP.^h);      %%%% \alpha

% Reaction rates
vC1 = kcat1 * T1 * C1 / (KT + C1);
vC2 = kcat2 * T2 * C2 / (KT + C2);
vM = kr * R * M / (Kr + M);
RPC1 = C1^2 / (gamma^2 + C1^2);
RP_RP = alpha^2 / (alpha^2 + RP^2);

% Differential equations
dC1_dt = -vC1;
dC2_dt = -vC2;
dM_dt = vC1 + vC2 - vM;
dQ_dt = (beta_p / (nQ)) * vM - kde * Q;
dR_dt = (beta_p / (nR)) * vM - kde * R;
dT1_dt = (beta_p / (nT1)) * vM - kdT1 * T1;
dT2_dt = (beta_p / (nT2)) * vM * RP_RP - kdT2 * T2;
dRP_dt = (beta_p / (nRP)) * vM * RPC1 - kdRP * RP;

% Pack derivatives into a column vector
dS = [dC1_dt; dC2_dt; dM_dt; dQ_dt; dR_dt; dT1_dt; dT2_dt; dRP_dt];

end
