%% Step 7: DMPC with communicated predecessor input
% Compared with Step 6, the current predecessor input u_1(k) is known by
% vehicle 2 before its local MPC optimization is solved.
% Thus x_(k+1|k) in X_R minus E*u_1(k), not X_R minus D.

clear; clc; close all;
yalmip('clear');

requiredFiles = {'terminal_set_step2.mat', 'robust_sets_step4.mat'};
for fileIndex = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{fileIndex})
        error('Missing %s. Run the earlier steps first.', requiredFiles{fileIndex});
    end
end

load('terminal_set_step2.mat', ...
    'A', 'B', 'P', 'Tset', ...
    'xmin', 'xmax', 'umin', 'umax', 'Qx', 'Qu', 'Ts');
load('robust_sets_step4.mat', 'E', 'Xrobust', 'uPrevMin', 'uPrevMax');

N = 11;
Nsim = 60;

Tset.minHRep();
Hterminal = Tset.A;
hterminal = Tset.b;

Xrobust.minHRep();
Hrobust = Xrobust.A;
hrobust = Xrobust.b;

%% The same predecessor input used in Step 6
uPredecessor = zeros(1, Nsim);
uPredecessor(1:8)   =  1.5;
uPredecessor(9:16)  = -2.0;
uPredecessor(17:24) =  0.8;

assert(all(uPredecessor >= uPrevMin) && all(uPredecessor <= uPrevMax), ...
    'The leader input must remain inside U^1.');

%% DMPC closed-loop simulation
x = zeros(2, Nsim+1);
uFollower = zeros(1, Nsim);
x(:,1) = [20; 5];

ops = sdpsettings('solver', 'quadprog', 'verbose', 0);

for k = 1:Nsim
    % Communication event: vehicle 2 receives the actual current u_1(k).
    uPredecessorKnown = uPredecessor(k);

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

    % Time-varying DMPC robustness constraint from Remark 3.6:
    % x_(k+1|k) + E*u_1(k) belongs to X_R.
    constraints = [constraints, ...
        Hrobust*(xPred(:,2) + E*uPredecessorKnown) <= hrobust, ...
        xmin <= xPred(:,N+1) <= xmax, ...
        Hterminal*xPred(:,N+1) <= hterminal];

    objective = objective + xPred(:,N+1)'*P*xPred(:,N+1);

    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem ~= 0
        error('DMPC became infeasible at k=%d: %s', k, diagnostics.info);
    end

    uFollower(k) = value(uPred(1));
    x(:,k+1) = A*x(:,k) + B*uFollower(k) + E*uPredecessorKnown;
end

fprintf('All DMPC optimization problems were feasible.\n');
fprintf('Final state: e_p = %.4g m, e_v = %.4g m/s\n', ...
    x(1,end), x(2,end));

%% DMPC behavior
figure('Color', 'w', 'Name', 'Step 7 - DMPC under communicated input');
tState = 0:Ts:Nsim*Ts;
tInput = 0:Ts:(Nsim-1)*Ts;

subplot(4,1,1)
plot(tState, x(1,:), 'LineWidth', 1.5); grid on
yline(xmin(1), '--r', 'Safety lower bound');
ylabel('e_p [m]')
title('DMPC with communicated predecessor input')

subplot(4,1,2)
plot(tState, x(2,:), 'LineWidth', 1.5); grid on
ylabel('e_v [m/s]')

subplot(4,1,3)
stairs(tInput, uPredecessor, 'LineWidth', 1.5); grid on
yline(uPrevMin, '--r'); yline(uPrevMax, '--r');
ylabel('u_1 [m/s^2]')
title('Communicated predecessor input')

subplot(4,1,4)
stairs(tInput, uFollower, 'LineWidth', 1.5); grid on
yline(umin, '--r'); yline(umax, '--r');
xlabel('time [s]')
ylabel('u_2 [m/s^2]')
title('Follower DMPC input')

%% Direct comparison with the RMPC simulation, if Step 6 has been run
if isfile('rmpc_simulation_step6.mat')
    rmpc = load('rmpc_simulation_step6.mat', 'x', 'uFollower', 'uPredecessor');

    if isequal(rmpc.uPredecessor, uPredecessor)
        figure('Color', 'w', 'Name', 'Step 7 - RMPC versus DMPC');
        subplot(3,1,1)
        plot(tState, rmpc.x(1,:), '--', 'LineWidth', 1.2); hold on
        plot(tState, x(1,:), 'LineWidth', 1.5); grid on
        ylabel('e_p [m]')
        legend('RMPC', 'DMPC', 'Location', 'best')
        title('Same plant and same predecessor input')

        subplot(3,1,2)
        plot(tState, rmpc.x(2,:), '--', 'LineWidth', 1.2); hold on
        plot(tState, x(2,:), 'LineWidth', 1.5); grid on
        ylabel('e_v [m/s]')
        legend('RMPC', 'DMPC', 'Location', 'best')

        subplot(3,1,3)
        stairs(tInput, rmpc.uFollower, '--', 'LineWidth', 1.2); hold on
        stairs(tInput, uFollower, 'LineWidth', 1.5); grid on
        xlabel('time [s]')
        ylabel('u_2 [m/s^2]')
        legend('RMPC', 'DMPC', 'Location', 'best')
    end
end

save('dmpc_simulation_step7.mat', 'x', 'uFollower', 'uPredecessor', ...
    'N', 'Nsim', 'Ts');
