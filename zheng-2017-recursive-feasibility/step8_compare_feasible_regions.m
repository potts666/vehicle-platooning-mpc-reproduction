%% Step 8: grid approximation of RMPC and DMPC feasible regions
% This reproduces the qualitative point of Figure 4 in the paper.  A grid
% point is marked feasible if the corresponding finite-horizon MPC problem
% has a solution.  The RMPC set is static; the DMPC set is shown for one
% communicated predecessor input u_1(k) = 1.5 m/s^2.

clear; clc; close all;
yalmip('clear');

requiredFiles = {'terminal_set_step2.mat', 'robust_sets_step4.mat', ...
                 'rmpc_constraint_step5.mat'};
for fileIndex = 1:numel(requiredFiles)
    if ~isfile(requiredFiles{fileIndex})
        error('Missing %s. Run the earlier steps first.', requiredFiles{fileIndex});
    end
end
load('terminal_set_step2.mat', ...
    'A', 'B', 'P', 'Tset', ...
    'xmin', 'xmax', 'umin', 'umax', 'Qx', 'Qu');
load('robust_sets_step4.mat', 'E', 'Xrobust', 'uPrevMin', 'uPrevMax');
load('rmpc_constraint_step5.mat', 'XRminusD');

N = 11;
uPredecessorKnown = 1.5;  % same value used during the first 8 s of Step 7

assert(uPredecessorKnown >= uPrevMin && uPredecessorKnown <= uPrevMax, ...
    'The selected communicated input must lie in U^1.');

Tset.minHRep();
Hterminal = Tset.A;
hterminal = Tset.b;

XRminusD.minHRep();
Hrm = XRminusD.A;
hrm = XRminusD.b;

Xrobust.minHRep();
Hdmpc = Xrobust.A;
hdmpc = Xrobust.b;

%% Build two parameterized feasibility problems once
rmpcController = buildFeasibilityController( ...
    A, B, E, Qx, Qu, P, N, xmin, xmax, umin, umax, ...
    Hterminal, hterminal, Hrm, hrm, 'rmpc', 0);

dmpcController = buildFeasibilityController( ...
    A, B, E, Qx, Qu, P, N, xmin, xmax, umin, umax, ...
    Hterminal, hterminal, Hdmpc, hdmpc, 'dmpc', uPredecessorKnown);

%% Scan the physical state-constraint rectangle
% Decrease the steps below if you want a faster but coarser figure.
epGrid = xmin(1):2:xmax(1);
evGrid = xmin(2):1:xmax(2);

feasibleRMPC = false(numel(evGrid), numel(epGrid));
feasibleDMPC = false(numel(evGrid), numel(epGrid));

fprintf('Scanning %d state-grid points for each controller...\n', ...
    numel(epGrid)*numel(evGrid));

for row = 1:numel(evGrid)
    for column = 1:numel(epGrid)
        x0 = [epGrid(column); evGrid(row)];

        [~, rmpcStatus] = rmpcController{x0};
        [~, dmpcStatus] = dmpcController{x0};

        feasibleRMPC(row,column) = (rmpcStatus == 0);
        feasibleDMPC(row,column) = (dmpcStatus == 0);
    end

    if mod(row,5) == 0 || row == numel(evGrid)
        fprintf('  completed %d of %d velocity-error rows\n', row, numel(evGrid));
    end
end

rmpcCount = nnz(feasibleRMPC);
dmpcCount = nnz(feasibleDMPC);
inclusionViolations = nnz(feasibleRMPC & ~feasibleDMPC);

fprintf('RMPC feasible grid points: %d\n', rmpcCount);
fprintf('DMPC feasible grid points: %d\n', dmpcCount);
fprintf('Grid points feasible for RMPC but not DMPC: %d\n', inclusionViolations);

%% Plot only the boundaries: RMPC dashed red, DMPC solid blue
figure('Color', 'w', 'Name', 'Step 8 - feasible-region comparison');
contour(epGrid, evGrid, double(feasibleRMPC), [0.5 0.5], ...
    '--r', 'LineWidth', 2);
hold on
contour(epGrid, evGrid, double(feasibleDMPC), [0.5 0.5], ...
    '-b', 'LineWidth', 2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
grid on
xlabel('e_p [m]')
ylabel('e_v [m/s]')
title(sprintf(['Finite-horizon feasible regions, ', ...
    'DMPC at communicated u_1 = %.1f m/s^2'], uPredecessorKnown))
legend('RMPC feasible-region boundary', 'DMPC feasible-region boundary', ...
    'Origin', 'Location', 'best')

save('feasible_regions_step8.mat', ...
    'epGrid', 'evGrid', 'feasibleRMPC', 'feasibleDMPC', ...
    'uPredecessorKnown', 'rmpcCount', 'dmpcCount', 'inclusionViolations');

%% Local function: build the finite-horizon MPC feasibility checker
function controller = buildFeasibilityController( ...
    A, B, E, Qx, Qu, P, N, xmin, xmax, umin, umax, ...
    Hterminal, hterminal, Hrobust, hrobust, mode, uPredecessorKnown)

    x0 = sdpvar(2,1);
    xPred = sdpvar(2,N+1,'full');
    uPred = sdpvar(1,N,'full');

    constraints = [xPred(:,1) == x0];
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

    if strcmp(mode, 'rmpc')
        % x_(k+1|k) belongs to X_R minus D.
        constraints = [constraints, Hrobust*xPred(:,2) <= hrobust];
    elseif strcmp(mode, 'dmpc')
        % x_(k+1|k) + E*u_1(k) belongs to X_R.
        constraints = [constraints, ...
            Hrobust*(xPred(:,2) + E*uPredecessorKnown) <= hrobust];
    else
        error('Unknown controller mode.');
    end

    constraints = [constraints, ...
        xmin <= xPred(:,N+1) <= xmax, ...
        Hterminal*xPred(:,N+1) <= hterminal];
    objective = objective + xPred(:,N+1)'*P*xPred(:,N+1);

    ops = sdpsettings('solver', 'quadprog', 'verbose', 0);
    controller = optimizer(constraints, objective, ops, x0, uPred(1));
end
