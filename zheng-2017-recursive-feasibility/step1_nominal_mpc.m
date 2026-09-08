%% Step 1: nominal MPC for one follower vehicle
% This is deliberately simpler than Problems 2.7 and 3.3 in the paper.
% It verifies the discrete model, LQR terminal cost, YALMIP, and QUADPROG
% before disturbance sets, RCI sets, RMPC, or DMPC are introduced.

clear; clc; close all;
yalmip('clear');

%% Parameters from the paper
Ts = 1;                     % sampling time [s]
h = 1;                      % velocity-dependent spacing policy
ds = 4;                     % standstill desired distance [m]

A = [1 Ts; 0 1];
B = [-h*Ts - Ts^2/2; -Ts];

xmin = [-ds; -15];          % [minimum position error; minimum speed error]
xmax = [120;  15];
umin = -5;                  % temporary physical input bound [m/s^2]
umax =  3;

Qx = eye(2);
Qu = 1;
N  = 11;                    % a reasonable starting horizon for this paper
Nsim = 60;

%% LQR terminal controller and terminal cost
% MATLAB uses u = -K*x. P defines the terminal quadratic cost x'*P*x.
[K, P, closedLoopPoles] = dlqr(A, B, Qx, Qu);
Acl = A - B*K;

assert(rank(ctrb(A,B)) == 2, 'The selected model is not controllable.');
assert(max(abs(eig(Acl))) < 1, 'The LQR closed loop is not stable.');

fprintf('LQR gain K = [%g  %g]\n', K(1), K(2));
fprintf('Closed-loop poles = [%g  %g]\n', ...
    closedLoopPoles(1), closedLoopPoles(2));

%% Closed-loop simulation
% x = [distance error; relative-speed error].
% The initial condition is chosen inside X and is only for this first test.
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

    constraints = [constraints, xmin <= xPred(:,N+1) <= xmax];
    objective = objective + xPred(:,N+1)'*P*xPred(:,N+1);

    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem ~= 0
        error('MPC became infeasible at k=%d: %s', k, diagnostics.info);
    end

    uApplied(k) = value(uPred(1));
    x(:,k+1) = A*x(:,k) + B*uApplied(k);
end

%% Basic checks and plots
fprintf('Final state: e_p = %.4g m, e_v = %.4g m/s\n', ...
    x(1,end), x(2,end));

figure('Color', 'w', 'Name', 'Step 1 - nominal MPC');
tState = 0:Ts:Nsim*Ts;
tInput = 0:Ts:(Nsim-1)*Ts;

subplot(3,1,1)
plot(tState, x(1,:), 'LineWidth', 1.5); grid on
yline(-ds, '--r', 'Safety lower bound');
ylabel('e_p [m]')
title('Nominal MPC: distance error')

subplot(3,1,2)
plot(tState, x(2,:), 'LineWidth', 1.5); grid on
ylabel('e_v [m/s]')
title('Nominal MPC: relative-speed error')

subplot(3,1,3)
stairs(tInput, uApplied, 'LineWidth', 1.5); grid on
yline(umin, '--r'); yline(umax, '--r');
xlabel('time [s]')
ylabel('u [m/s^2]')
title('Applied acceleration')
