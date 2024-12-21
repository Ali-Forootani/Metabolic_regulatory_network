% Parameters


n_Q = 300;
n_RP = 300;
n_T1 = 400;
n_T2 = 1500;
n_R = 7459;
w = 100;

% Initial conditions
x0 = [500, 1000, 20, 0.1, 0.01, 0.01, 0.01, 0.01];  % C1, C2, M, Q, R, T1, T2, RP

% Time span
tspan = [0, 70];
t_eval = linspace(0, 70, 1000);  % Higher resolution for smooth curves
% Solving the ODE
[t, sol] = ode45(@myFunc_new, t_eval, x0);

% Calculating biomass
% Biomass = w * sol(:, 3) + w * n_Q * sol(:, 4) + w * n_R * sol(:, 5) + ...
%           w * n_T1 * sol(:, 6) + w * n_T2 * sol(:, 7) + w * n_RP * sol(:, 8);

% Plotting Concentration over Time
 figure;
% yyaxis left
 plot(t, sol(:, 1), 'g', 'LineWidth', 2);
% hold on;
% plot(t, sol(:, 2), 'b', 'LineWidth', 2);
% plot(t, sol(:, 3), 'r', 'LineWidth', 2);
% xlabel('Time (min)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Metabolite Concentration (mmol)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Concentration over Time', 'FontSize', 14, 'FontWeight', 'bold');
% legend('C_1', 'C_2', 'M');
% grid on;

% Plotting Biomass over Time
% figure;
% plot(t, Biomass, 'k--', 'LineWidth', 2);
% xlabel('Time (min)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Biomass (mg)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Biomass over Time', 'FontSize', 14, 'FontWeight', 'bold');
% grid on;

% Plotting Macromolecular Mass over Time
% figure;
% plot(t, w * n_Q * sol(:, 4), 'b', 'LineWidth', 2);
% hold on;
% plot(t, w * n_R * sol(:, 5), 'g', 'LineWidth', 2);
% plot(t, w * n_T1 * sol(:, 6), 'r', 'LineWidth', 2);
% plot(t, w * n_T2 * sol(:, 7), 'k', 'LineWidth', 2);
% plot(t, w * n_RP * sol(:, 8), 'c', 'LineWidth', 2);
% xlabel('Time (min)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Macromolecular Mass (mg)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Macromolecular Mass over Time', 'FontSize', 14, 'FontWeight', 'bold');
% legend('Q', 'R', 'T1', 'T2', 'RP');
% grid on;
% ODE function
function eqn = myFunc_new(t, s)
    C1 = s(1); C2 = s(2); M = s(3); Q = s(4);
    R = s(5); T1 = s(6); T2 = s(7); RP = s(8);
n_Q = 300;
n_RP = 300;
n_T1 = 400;
n_T2 = 1500;
n_R = 7459;
w = 100;
    % Parameters
    k_1 = 3000;  % kcatc1
    k_2 = 2000;  % kcatc2
    k_r = 1260;  % kcatQ, kcatR, kcatT1, kcatT2, kcatRP

    % Degradation constants
    k_de = 0.01;  % kdegQ, kdegR
    k_dT1 = 0.05;  % kdegT1
    k_dT2 = 0.05;  % kdegT2
    k_dRP = 0.2;  % kdegRP

    % Michaelis constants
    K_T = 1000;  % kMc1, kMc2
    K_r = 7;     % kMQ, kMR, kMT1, kMT2, kMRP

    % Hill function parameters
    h = 2;
    gamma = 20;
    alpha = 0.03;
    RPC1 = C1^h / (gamma^h + C1^h);
    RP_RP = alpha^h / (alpha^h + RP^h);

    % Reaction rates
    vc1 = k_1 * C1 * T1 / (K_T + C1);
    vc2 = k_2 * C2 * T2 / (K_T + C2);
    vM = k_r * M * R / (K_r + M);

    % ODEs
    eqn = zeros(8, 1);
    eqn(1) = -vc1;  % C1
    eqn(2) = -vc2;  % C2
    eqn(3) = vc1 + vc2 - vM;  % M
    eqn(4) = (1 / n_Q) * vM - k_de * Q;  % Q
    eqn(5) = (1 /  n_R) * vM - k_de * R;  % R
    eqn(6) = (1 /  n_T1) * vM - k_dT1 * T1;  % T1
    eqn(7) = (1 /  n_T2) * vM * RP_RP - k_dT2 * T2;  % T2
    eqn(8) = (1 / n_RP) * vM * RPC1 - k_dRP * RP;  % RP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coffeicient

function eqn = sindy(t, s)
    C1 = s(1); C2 = s(2); M = s(3); Q = s(4);
    R = s(5); T1 = s(6); T2 = s(7); RP = s(8);
n_Q = 300;
n_RP = 3000;
n_T1 = 400;
n_T2 = 1500;
n_R = 7459;
w = 100;
    % Parameters
    k_1 = 3000;  % kcatc1
    k_2 = 2000;  % kcatc2
    k_r = 1260;  % kcatQ, kcatR, kcatT1, kcatT2, kcatRP

    % Degradation constants
    k_de = 0.01;  % kdegQ, kdegR
    k_dT1 = 0.05;  % kdegT1
    k_dT2 = 0.05;  % kdegT2
    k_dRP = 0.2;  % kdegRP

    % Michaelis constants
    K_T = 1000;  % kMc1, kMc2
    K_r = 7;     % kMQ, kMR, kMT1, kMT2, kMRP

    % Hill function parameters
    h = 2;
    gamma = 20;
    alpha = 0.03;
    RPC1 = C1^h / (gamma^h + C1^h);
    RP_RP = alpha^h / (alpha^h + RP^h);

    % Reaction rates
    vc1 = k_1 * C1 * T1 / (K_T + C1);
    vc2 = k_2 * C2 * T2 / (K_T + C2);
    vM = k_r * M * R / (K_r + M);

    % ODEs
    eqn = zeros(8, 1);
    eqn(1) = -vc1;  % C1
    eqn(2) = -vc2;  % C2
    eqn(3) = vc1 + vc2 - vM;  % M
    eqn(4) = (1 / n_Q) * vM - k_de * Q;  % Q
    eqn(5) = (1 /  n_R) * vM - k_de * R;  % R
    eqn(6) = (1 /  n_T1) * vM - k_dT1 * T1;  % T1
    eqn(7) = (1 /  n_T2) * vM * RP_RP - k_dT2 * T2;  % T2
    eqn(8) = (1 / n_RP) * vM * RPC1 - k_dRP * RP;  % RP
end