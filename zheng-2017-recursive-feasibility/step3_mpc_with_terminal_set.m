%% Step 3: nominal MPC with the terminal-set constraint x_N|k in T^i
% This is the nominal part of the paper's local MPC formulation.  It still
% does not include the predecessor-input disturbance or robust constraint.

clear; clc; close all;
yalmip('clear');

if ~isfile('terminal_set_step2.mat')
    error(['terminal_set_step2.mat was not found. ', ...
        'Run step2_terminal_set.m first.']);
end

load('terminal_set_step2.mat', ...
    'A', 'B', 'K', 'P', 'Tset', ...
    'xmin', 'xmax', 'umin', 'umax', 'Qx', 'Qu', 'Ts');

N = 11;
Nsim = 60;

% Tset = {x | Hterminal*x <= hterminal}
Tset.minHRep();
Hterminal = Tset.A;
hterminal = Tset.b;

%% Closed-loop simulation
x = zeros(2, Nsim+1);
uApplied = zeros(1, Nsim);
x(:,1) = [20; 5];

ops = sdpsettings('solver', 'quadprog', 'verbose', 0);

for k = 1:Nsim
    xPred = sdpvar(2, N+1, 'full');
    uPred = sdpvar(1, N, 'full');

    constraints = [xPred(:,1) == x(:,k)];
    objective = 0;

    for j = 1:N
        constraints = [constraints, ...
            xPred(:,j+1) == A*xPred(:,j) + B*uPred(j), ...
            xmin <= xPred(:,j) <= xmax, ...
            umin <= uPred(j) <= umax];

        objective = objective ...
            + xPred(:,j)'*Qx*xPred(:,j) ...
            + uPred(j)'*Qu*uPred(j);
    end

    constraints = [constraints, ...
        xmin <= xPred(:,N+1) <= xmax, ...
        Hterminal*xPred(:,N+1) <= hterminal];

    objective = objective + xPred(:,N+1)'*P*xPred(:,N+1);

    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem ~= 0
        error('MPC became infeasible at k=%d: %s', k, diagnostics.info);
    end

    uApplied(k) = value(uPred(1));
    x(:,k+1) = A*x(:,k) + B*uApplied(k);
end

fprintf('Final state: e_p = %.4g m, e_v = %.4g m/s\n', ...
    x(1,end), x(2,end));
fprintf('All terminal-set constrained MPC problems were feasible.\n');

%% Show the closed-loop path approaching the terminal set
figure('Color', 'w', 'Name', 'Step 3 - MPC with terminal set');
Tset.plot('color', [0 0.45 0.74], 'linewidth', 2);
hold on
plot(x(1,:), x(2,:), '-o', 'Color', [0.85 0.33 0.10], ...
    'MarkerSize', 3, 'LineWidth', 1.2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
grid on
xlabel('e_p [m]')
ylabel('e_v [m/s]')
title('Closed-loop MPC trajectory and terminal set T^i')
legend('Terminal set T^i', 'Closed-loop state trajectory', 'Origin', ...
    'Location', 'best')

figure('Color', 'w', 'Name', 'Step 3 - terminal-set constrained MPC');
tState = 0:Ts:Nsim*Ts;
tInput = 0:Ts:(Nsim-1)*Ts;

subplot(3,1,1)
plot(tState, x(1,:), 'LineWidth', 1.5); grid on
yline(xmin(1), '--r', 'Safety lower bound');
ylabel('e_p [m]')

subplot(3,1,2)
plot(tState, x(2,:), 'LineWidth', 1.5); grid on
ylabel('e_v [m/s]')

subplot(3,1,3)
stairs(tInput, uApplied, 'LineWidth', 1.5); grid on
yline(umin, '--r'); yline(umax, '--r');
xlabel('time [s]')
ylabel('u [m/s^2]')
