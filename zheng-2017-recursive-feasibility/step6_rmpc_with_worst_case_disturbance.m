%% Step 6: RMPC with the static robustness constraint X_R^i minus D^i
% The optimizer uses the nominal prediction model x+ = A*x+B*u.
% The simulated plant uses x+ = A*x+B*u+E*u_predecessor.
% The static first-step constraint protects against every allowed u_predecessor.

clear; clc; close all;
yalmip('clear');

requiredFiles = {'terminal_set_step2.mat', 'rmpc_constraint_step5.mat', ...
                 'robust_sets_step4.mat'};
for fileIndex = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{fileIndex})
        error('Missing %s. Run the earlier steps first.', requiredFiles{fileIndex});
    end
end

load('terminal_set_step2.mat', ...
    'A', 'B', 'P', 'Tset', ...
    'xmin', 'xmax', 'umin', 'umax', 'Qx', 'Qu', 'Ts');
load('robust_sets_step4.mat', 'E', 'uPrevMin', 'uPrevMax');
load('rmpc_constraint_step5.mat', 'XRminusD');

N = 11;
Nsim = 60;

Tset.minHRep();
Hterminal = Tset.A;
hterminal = Tset.b;

XRminusD.minHRep();
Hrobust = XRminusD.A;
hrobust = XRminusD.b;

%% A deterministic leader input that stays within U^1
% RMPC does NOT receive this signal while optimizing.  It only knows
% that every value lies in [uPrevMin, uPrevMax].
uPredecessor = zeros(1, Nsim);
uPredecessor(1:8)   =  1.5;
uPredecessor(9:16)  = -2.0;
uPredecessor(17:24) =  0.8;

assert(all(uPredecessor >= uPrevMin) && all(uPredecessor <= uPrevMax), ...
    'The leader input must remain inside U^1.');

%% RMPC closed-loop simulation
x = zeros(2, Nsim+1);
uFollower = zeros(1, Nsim);
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

    % Paper's RMPC robustness constraint:
    % x_(k+1|k) must lie in X_R minus D.
    constraints = [constraints, ...
        Hrobust*xPred(:,2) <= hrobust, ...
        xmin <= xPred(:,N+1) <= xmax, ...
        Hterminal*xPred(:,N+1) <= hterminal];

    objective = objective + xPred(:,N+1)'*P*xPred(:,N+1);

    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem ~= 0
        error('RMPC became infeasible at k=%d: %s', k, diagnostics.info);
    end

    uFollower(k) = value(uPred(1));

    % Actual plant: predecessor input enters as an unknown disturbance.
    x(:,k+1) = A*x(:,k) + B*uFollower(k) + E*uPredecessor(k);
end

fprintf('All RMPC optimization problems were feasible.\n');
fprintf('Final state: e_p = %.4g m, e_v = %.4g m/s\n', ...
    x(1,end), x(2,end));

%% Plot the actual closed-loop behavior
figure('Color', 'w', 'Name', 'Step 6 - RMPC under predecessor input');
tState = 0:Ts:Nsim*Ts;
tInput = 0:Ts:(Nsim-1)*Ts;

subplot(4,1,1)
plot(tState, x(1,:), 'LineWidth', 1.5); grid on
yline(xmin(1), '--r', 'Safety lower bound');
ylabel('e_p [m]')
title('RMPC under an unknown-but-bounded predecessor input')

subplot(4,1,2)
plot(tState, x(2,:), 'LineWidth', 1.5); grid on
ylabel('e_v [m/s]')

subplot(4,1,3)
stairs(tInput, uPredecessor, 'LineWidth', 1.5); grid on
yline(uPrevMin, '--r'); yline(uPrevMax, '--r');
ylabel('u_1 [m/s^2]')
title('Actual predecessor input (not communicated to RMPC)')

subplot(4,1,4)
stairs(tInput, uFollower, 'LineWidth', 1.5); grid on
yline(umin, '--r'); yline(umax, '--r');
xlabel('time [s]')
ylabel('u_2 [m/s^2]')
title('Follower RMPC input')

save('rmpc_simulation_step6.mat', 'x', 'uFollower', 'uPredecessor', ...
    'N', 'Nsim', 'Ts');
