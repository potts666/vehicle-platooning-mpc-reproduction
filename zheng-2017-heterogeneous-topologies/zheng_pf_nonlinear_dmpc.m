%% Zheng et al. (2017): nonlinear DMPC platoon simulation
% One leader and seven heterogeneous followers.  This script implements
% (1), (13)--(18) for PF, PLF, TPF, or TPLF communication topologies.
%
% Requirements: MATLAB Optimization Toolbox (fmincon).
% The paper specifies m_i, tau_i, C_A,i, R_i and the cost weights.  It does
% not list eta_i or f_i in Table I; the values below are explicit modelling
% choices for this reproducible simulation.

clearvars -except TsimOverride; clc; close all;

if exist('fmincon', 'file') ~= 2
    error('This script requires fmincon from MATLAB Optimization Toolbox.');
end

%% Paper simulation setting
dt = 0.1;                 % Paper: Delta t = 0.1 s
Np = 20;                  % Paper: predictive horizon
Nf = 7;                   % Seven followers, plus one leader
Tsim = 6;                 % Fig. 4 displays 0--6 s
if exist('TsimOverride', 'var')
    Tsim = TsimOverride;  % Optional short-run check, e.g. TsimOverride = 0.5
end
Nsim = round(Tsim / dt);
d0 = 20;                  % Constant desired spacing [m]
g = 9.81;
v0 = 20;                  % Initial leader/follower speed [m/s]
aMax = 6;                 % Paper acceleration limit [m/s^2]
aMin = -6;                % Paper deceleration limit [m/s^2]

% Change only this line to select a topology:
% 'PF'   = predecessor following
% 'PLF'  = predecessor-leader following
% 'TPF'  = two-predecessor following
% 'TPLF' = two-predecessor-leader following
topology = 'TPLF';
[predecessors, leaderPinned] = topology_definition(topology, Nf);

% Table I of the paper: [mass, powertrain lag, aerodynamic coefficient, R].
tableI = [ ...
    1035.7, 0.51, 0.99, 0.30; ...
    1849.1, 0.75, 1.15, 0.38; ...
    1934.0, 0.78, 1.17, 0.39; ...
    1678.7, 0.70, 1.12, 0.37; ...
    1757.7, 0.73, 1.13, 0.38; ...
    1743.1, 0.72, 1.13, 0.37; ...
    1392.2, 0.62, 1.06, 0.34];

% Explicit choices not tabulated in the paper.
eta = 0.90 * ones(Nf, 1); % Driveline mechanical efficiency
fRoll = 0.015 * ones(Nf, 1); % Rolling-resistance coefficient

params = repmat(struct('m', [], 'tau', [], 'CA', [], 'R', [], ...
    'eta', [], 'fRoll', [], 'umin', [], 'umax', []), Nf, 1);
for i = 1:Nf
    params(i).m = tableI(i, 1);
    params(i).tau = tableI(i, 2);
    params(i).CA = tableI(i, 3);
    params(i).R = tableI(i, 4);
    params(i).eta = eta(i);
    params(i).fRoll = fRoll(i);

    % Convert the +/-6 m/s^2 paper acceleration range into fixed command-
    % torque bounds around 22 m/s.  The actual vehicle still has the
    % first-order torque lag in (1).
    force22 = params(i).CA * 22^2 + params(i).m * g * params(i).fRoll;
    params(i).umin = params(i).R / params(i).eta * ...
        (params(i).m * aMin + force22);
    params(i).umax = params(i).R / params(i).eta * ...
        (params(i).m * aMax + force22);
end

%% Table II weights. F_i = 10 I_2, G_i = 5 I_2 for every vehicle that
% has predecessor information, R_i = 1, and Q_i is nonzero exactly for
% leader-pinned vehicles.
F = 10 * eye(2);
G = 5 * eye(2);
Q = zeros(2, 2, Nf);
for i = 1:Nf
    if leaderPinned(i)
        Q(:, :, i) = 10 * eye(2);
    end
end
Rweight = 1;

%% Initial state: the desired formation at 20 m/s
% x_i = [position; speed; actual torque].
x = zeros(3, Nf);
for i = 1:Nf
    x(:, i) = [-i * d0; v0; equilibrium_torque(v0, params(i), g)];
end

stateHistory = zeros(3, Nf, Nsim + 1);
stateHistory(:, :, 1) = x;
uApplied = zeros(Nf, Nsim);
leaderHistory = zeros(2, Nsim + 1);
leaderHistory(:, 1) = leader_output(0);

% Initialization in (15): each follower assumes constant-speed motion.
uAssumed = zeros(Nf, Np);
yAssumed = zeros(Nf, 2, Np + 1);
for i = 1:Nf
    uAssumed(i, :) = equilibrium_torque(v0, params(i), g);
    xInit = rollout(x(:, i), uAssumed(i, :).', params(i), dt, g);
    yAssumed(i, :, :) = reshape(xInit(1:2, :), [1, 2, Np + 1]);
end

solverOptions = optimoptions('fmincon', ...
    'Algorithm', 'sqp', 'Display', 'off', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 30000, ...
    'ConstraintTolerance', 1e-7, 'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-9, 'ScaleProblem', 'obj-and-constr');
% Terminal position, speed, and torque have different physical units and
% magnitudes.  The nonlinear constraints are therefore nondimensionalized
% before they are passed to fmincon.  This does not change their feasible set.
feasibilityAcceptanceTol = 5e-6;

%% Iteration of DMPC, equations (16)--(18)
for tIndex = 1:Nsim
    tNow = (tIndex - 1) * dt;

    % All local problems must use the same old assumed trajectories.
    yAssumedOld = yAssumed;
    uStar = zeros(Nf, Np);
    xStar = cell(Nf, 1);

    % Future leader output is known only to the pinned first follower.
    leaderFuture = zeros(2, Np + 1);
    for k = 0:Np
        leaderFuture(:, k + 1) = leader_output(tNow + k * dt);
    end

    for i = 1:Nf
        ySelfAssumed = squeeze(yAssumedOld(i, :, :));

        % Neighbor outputs used in the G_i term of (14).  A vehicle i is
        % (i-j)*d0 behind predecessor j.  The leader is not included here:
        % its direct information is represented by Q_i when i is pinned.
        neighborOutputs = cell(numel(predecessors{i}), 1);
        neighborOffsets = cell(numel(predecessors{i}), 1);
        for r = 1:numel(predecessors{i})
            j = predecessors{i}(r);
            neighborOutputs{r} = squeeze(yAssumedOld(j, :, :));
            neighborOffsets{r} = [-(i - j) * d0; 0];
        end

        % Terminal constraint (13d): average the offset-corrected terminal
        % assumed outputs of all members of I_i = N_i union P_i.
        terminalReferences = {};
        if leaderPinned(i)
            terminalReferences{end + 1} = leaderFuture(:, end) + [-i * d0; 0];
        end
        for r = 1:numel(predecessors{i})
            terminalReferences{end + 1} = ...
                neighborOutputs{r}(:, end) + neighborOffsets{r};
        end
        terminalTarget = mean(cat(2, terminalReferences{:}), 2);

        % Only vehicle 1 uses a direct desired-trajectory penalty.  For the
        % remaining PF vehicles Q_i = 0, exactly as in Table II.
        yDesired = leaderFuture(:, 1:Np) + [-i * d0; 0];
        terminalScales = [d0; 22; max(abs([params(i).umin, params(i).umax]))];

        u0 = uAssumed(i, :).';
        objective = @(u) local_cost(u, x(:, i), params(i), dt, g, ...
            ySelfAssumed, neighborOutputs, neighborOffsets, yDesired, ...
            Q(:, :, i), F, G, Rweight);
        equalityConstraints = @(u) terminal_constraints(u, x(:, i), ...
            params(i), dt, g, terminalTarget, terminalScales);

        [uStar(i, :), ~, exitflag, output] = fmincon(objective, u0, ...
            [], [], [], [], params(i).umin * ones(Np, 1), ...
            params(i).umax * ones(Np, 1), equalityConstraints, solverOptions);
        [cCheck, ceqCheck, rawCeqCheck] = terminal_constraints( ...
            uStar(i, :).', x(:, i), params(i), dt, g, terminalTarget, terminalScales);
        maxViolation = max([0; cCheck(:); abs(ceqCheck(:))]);
        if maxViolation > feasibilityAcceptanceTol
            error(['Local MPC for vehicle %d failed at t = %.2f s. ', ...
                'Maximum scaled terminal residual = %.3e; ', ...
                'raw residuals [position, speed, torque] = [%.3e, %.3e, %.3e]. ', ...
                'fmincon message: %s'], i, tNow, maxViolation, ...
                rawCeqCheck(1), rawCeqCheck(2), rawCeqCheck(3), output.message);
        elseif exitflag <= 0
            warning(['Vehicle %d at t = %.2f s returned a near-feasible ', ...
                'point (residual %.3e): %s'], ...
                i, tNow, maxViolation, output.message);
        end
        xStar{i} = rollout(x(:, i), uStar(i, :).', params(i), dt, g);
    end

    % Equations (17)--(18): shift optimal inputs, append equilibrium torque,
    % and generate the assumed trajectory to communicate at t+1.
    uAssumedNext = zeros(Nf, Np);
    yAssumedNext = zeros(Nf, 2, Np + 1);
    for i = 1:Nf
        uAssumedNext(i, 1:Np-1) = uStar(i, 2:Np);
        uAssumedNext(i, Np) = equilibrium_torque(xStar{i}(2, end), params(i), g);
        xAssumedNext = rollout(xStar{i}(:, 2), uAssumedNext(i, :).', ...
            params(i), dt, g);
        yAssumedNext(i, :, :) = reshape(xAssumedNext(1:2, :), [1, 2, Np + 1]);
    end

    % Apply only the first optimal input to the physical vehicles.
    for i = 1:Nf
        uApplied(i, tIndex) = uStar(i, 1);
        x(:, i) = vehicle_step(x(:, i), uApplied(i, tIndex), params(i), dt, g);
    end
    stateHistory(:, :, tIndex + 1) = x;
    leaderHistory(:, tIndex + 1) = leader_output(tNow + dt);
    uAssumed = uAssumedNext;
    yAssumed = yAssumedNext;
end

%% Spacing-error verification and figures
time = 0:dt:Tsim;
spacingError = zeros(Nf, Nsim + 1);
for k = 1:Nsim+1
    for i = 1:Nf
        predecessorPosition = leaderHistory(1, k);
        if i > 1
            predecessorPosition = stateHistory(1, i - 1, k);
        end
        spacingError(i, k) = predecessorPosition - stateHistory(1, i, k) - d0;
    end
end

fprintf('%s nonlinear DMPC simulation completed.\n', topology);
fprintf('Maximum absolute spacing error: %.3f m\n', max(abs(spacingError), [], 'all'));
fprintf('Minimum actual adjacent gap: %.3f m\n', d0 + min(spacingError, [], 'all'));

figure('Color', 'w', 'Name', sprintf('Zheng 2017 %s nonlinear DMPC', topology));
subplot(3, 1, 1)
plot(time, spacingError, 'LineWidth', 1.25); grid on
yline(0, 'k:'); ylabel('spacing error [m]')
title(sprintf('%s topology: spacing errors', topology))
legend(compose('vehicle %d', 1:Nf), 'Location', 'eastoutside')

subplot(3, 1, 2)
plot(time, leaderHistory(2, :), 'k--', 'LineWidth', 1.8); hold on
for i = 1:Nf
    plot(time, squeeze(stateHistory(2, i, :)), 'LineWidth', 1.1);
end
grid on; ylabel('speed [m/s]')
title('Leader and follower speeds')
legend(["leader", compose('vehicle %d', 1:Nf)], 'Location', 'eastoutside')

subplot(3, 1, 3)
stairs(time(1:end-1), uApplied.', 'LineWidth', 1.1); grid on
xlabel('time [s]'); ylabel('command torque [N m]')
title('Applied first element of each optimal input sequence')

%% Local functions
function xNext = vehicle_step(x, u, p, dt, g)
    dragAndRolling = p.CA * x(2)^2 + p.m * g * p.fRoll;
    xNext = [x(1) + dt * x(2); ...
             x(2) + dt / p.m * (p.eta / p.R * x(3) - dragAndRolling); ...
             x(3) + dt / p.tau * (u - x(3))];
end

function X = rollout(x0, u, p, dt, g)
    NpLocal = numel(u);
    X = zeros(3, NpLocal + 1);
    X(:, 1) = x0;
    for k = 1:NpLocal
        X(:, k + 1) = vehicle_step(X(:, k), u(k), p, dt, g);
    end
end

function cost = local_cost(u, x0, p, dt, g, ySelfA, neighborOutputs, ...
        neighborOffsets, yDesired, Q, F, G, Rweight)
    X = rollout(x0, u, p, dt, g);
    cost = 0;
    for k = 1:numel(u)
        y = X(1:2, k);
        uEquilibrium = equilibrium_torque(X(2, k), p, g);
        eDesired = y - yDesired(:, k);
        eSelf = y - ySelfA(:, k);
        cost = cost + eDesired.' * Q * eDesired + ...
            Rweight * (u(k) - uEquilibrium)^2 + ...
            eSelf.' * F * eSelf;

        for r = 1:numel(neighborOutputs)
            eNeighbor = y - neighborOutputs{r}(:, k) - neighborOffsets{r};
            cost = cost + eNeighbor.' * G * eNeighbor;
        end
    end
end

function [c, ceq, rawCeq] = terminal_constraints(u, x0, p, dt, g, ...
        terminalTarget, terminalScales)
    X = rollout(x0, u, p, dt, g);
    yTerminal = X(1:2, end);
    c = [];
    rawCeq = [yTerminal - terminalTarget; ...
               X(3, end) - equilibrium_torque(X(2, end), p, g)];
    ceq = rawCeq ./ terminalScales;
end

function torque = equilibrium_torque(v, p, g)
    torque = p.R / p.eta * (p.CA * v^2 + p.m * g * p.fRoll);
end

function y = leader_output(t)
    % Continuous 20 -> 22 m/s ramp.  The displayed paper expression omits
    % '(t-1)' in its middle branch; this continuous form matches its text.
    if t <= 1
        v = 20;
        s = 20 * t;
    elseif t <= 2
        q = t - 1;
        v = 20 + 2 * q;
        s = 20 + 20 * q + q^2;
    else
        v = 22;
        s = 41 + 22 * (t - 2);
    end
    y = [s; v];
end

function [predecessors, leaderPinned] = topology_definition(topology, Nf)
    predecessors = cell(Nf, 1);
    leaderPinned = false(Nf, 1);

    switch upper(topology)
        case 'PF'
            % 0 -> 1 -> 2 -> ... -> Nf
            leaderPinned(1) = true;
            for i = 2:Nf
                predecessors{i} = i - 1;
            end

        case 'PLF'
            % Each vehicle follows its predecessor and knows the leader.
            leaderPinned(:) = true;
            for i = 2:Nf
                predecessors{i} = i - 1;
            end

        case 'TPF'
            % Each vehicle follows its two closest predecessors when present.
            leaderPinned(1:2) = true;
            predecessors{2} = 1;
            for i = 3:Nf
                predecessors{i} = [i - 1, i - 2];
            end

        case 'TPLF'
            % Two-predecessor following and direct leader information for all.
            leaderPinned(:) = true;
            predecessors{2} = 1;
            for i = 3:Nf
                predecessors{i} = [i - 1, i - 2];
            end

        otherwise
            error('Unknown topology: %s', topology);
    end
end
