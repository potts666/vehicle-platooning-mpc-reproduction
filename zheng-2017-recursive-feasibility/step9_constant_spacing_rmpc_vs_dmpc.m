%% Step 9: constant-spacing policy (h = 0), RMPC versus DMPC
% This is the qualitative counterpart of Section 4.1 / Figure 2 of the
% paper.  With C = diag(0.9,0.9), the static RMPC set X_R minus D may not
% contain the origin, while DMPC can use the current communicated input.

clear; clc; close all;
yalmip('clear');
mpt_init;

%% Constant-spacing model and common parameters
Ts = 1;
h = 0;
ds = 4;

A = [1 Ts; 0 1];
B = [-h*Ts - Ts^2/2; -Ts];
E = [Ts^2/2; Ts];

xmin = [-ds; -15];
xmax = [120;  15];
umin = -5;
umax =  3;
Qx = eye(2);
Qu = 1;
N = 11;
Nsim = 60;

C = diag([0.9 0.9]);
Upredecessor = [umin umax]*C;
uPrevMin = Upredecessor(1);
uPrevMax = Upredecessor(2);

X = Polyhedron('lb', xmin, 'ub', xmax);

%% Terminal ingredients and robust control invariant set
[K, P] = dlqr(A, B, Qx, Qu);
Tset = computeTerminalSet(A, B, K, xmin, xmax, umin, umax);
Xrobust = computeRCI(X, A, B, E, umin, umax, uPrevMin, uPrevMax);
XRminusD = pontryaginDifferenceWithInputSet( ...
    Xrobust, E, uPrevMin, uPrevMax);

origin = [0; 0];
originInRMPCConstraint = XRminusD.contains(origin);
originInteriorRMPCConstraint = all(XRminusD.b - XRminusD.A*origin > 1e-8);

fprintf('Constant-spacing model: B = [%g, %g]^T, E = [%g, %g]^T\n', ...
    B(1), B(2), E(1), E(2));
fprintf('Origin belongs to X_R minus D: %s\n', logicalText(originInRMPCConstraint));
fprintf('Origin strictly inside X_R minus D: %s\n', ...
    logicalText(originInteriorRMPCConstraint));

%% Same bounded predecessor input for both controllers
uPredecessor = zeros(1,Nsim);
uPredecessor(1:8)   =  1.5;
uPredecessor(9:16)  = -2.0;
uPredecessor(17:24) =  0.8;

assert(all(uPredecessor >= uPrevMin) && all(uPredecessor <= uPrevMax), ...
    'The predecessor input must remain within U^1.');

%% Simulate RMPC and DMPC from exactly the same initial condition
x0 = [20; 5];
[xRMPC, uRMPC] = simulateController( ...
    'rmpc', x0, uPredecessor, A, B, E, Qx, Qu, P, Tset, XRminusD, ...
    xmin, xmax, umin, umax, N);

[xDMPC, uDMPC] = simulateController( ...
    'dmpc', x0, uPredecessor, A, B, E, Qx, Qu, P, Tset, Xrobust, ...
    xmin, xmax, umin, umax, N);

fprintf('RMPC final state: e_p = %.4g m, e_v = %.4g m/s\n', ...
    xRMPC(1,end), xRMPC(2,end));
fprintf('DMPC final state: e_p = %.4g m, e_v = %.4g m/s\n', ...
    xDMPC(1,end), xDMPC(2,end));

%% Main comparison: Figure-2-like distance-error plot
tState = 0:Ts:Nsim*Ts;
tInput = 0:Ts:(Nsim-1)*Ts;

figure('Color', 'w', 'Name', 'Step 9 - constant spacing comparison');
plot(tState, xRMPC(1,:), '--', 'LineWidth', 1.8); hold on
plot(tState, xDMPC(1,:), 'LineWidth', 1.8);
yline(-ds, '--r', 'Safety lower bound');
grid on
xlabel('time [s]')
ylabel('e_p^2 [m]')
title('Constant spacing: RMPC versus DMPC')
legend('RMPC', 'DMPC', 'Location', 'best')

figure('Color', 'w', 'Name', 'Step 9 - constant spacing details');
subplot(3,1,1)
plot(tState, xRMPC(2,:), '--', 'LineWidth', 1.4); hold on
plot(tState, xDMPC(2,:), 'LineWidth', 1.4); grid on
ylabel('e_v^2 [m/s]')
legend('RMPC', 'DMPC', 'Location', 'best')

subplot(3,1,2)
stairs(tInput, uPredecessor, 'LineWidth', 1.4); grid on
yline(uPrevMin, '--r'); yline(uPrevMax, '--r');
ylabel('u_1 [m/s^2]')
title('Predecessor input')

subplot(3,1,3)
stairs(tInput, uRMPC, '--', 'LineWidth', 1.4); hold on
stairs(tInput, uDMPC, 'LineWidth', 1.4); grid on
yline(umin, '--r'); yline(umax, '--r');
xlabel('time [s]')
ylabel('u_2 [m/s^2]')
legend('RMPC', 'DMPC', 'Location', 'best')

save('constant_spacing_step9.mat', ...
    'A', 'B', 'E', 'K', 'P', 'Tset', 'Xrobust', 'XRminusD', ...
    'xRMPC', 'uRMPC', 'xDMPC', 'uDMPC', 'uPredecessor', ...
    'originInRMPCConstraint', 'originInteriorRMPCConstraint');

%% Local functions
function Tset = computeTerminalSet(A, B, K, xmin, xmax, umin, umax)
    Acl = A - B*K;
    H = [eye(2); -eye(2); -K; K];
    g = [xmax; -xmin; umax; -umin];
    Tset = Polyhedron('A', H, 'b', g);
    Tset.minHRep();

    for j = 1:100
        Tnext = Tset.intersect(Tset.invAffineMap(Acl));
        Tnext.minHRep();
        if Tnext == Tset
            Tset = Tnext;
            return
        end
        Tset = Tnext;
    end
    error('Terminal set did not converge.');
end
function Xrobust = computeRCI(X, A, B, E, umin, umax, uPrevMin, uPrevMax)
    Xrobust = X;
    for j = 1:100
        Pre = robustPredecessor(Xrobust, A, B, E, ...
            umin, umax, uPrevMin, uPrevMax);
        Xnext = X.intersect(Pre);
        Xnext.minHRep();
        if Xnext.isEmptySet()
            error('The robust control invariant set is empty.');
        end
        if Xnext == Xrobust
            Xrobust = Xnext;
            return
        end
        Xrobust = Xnext;
    end
    error('Robust control invariant set did not converge.');
end

function Pre = robustPredecessor(S, A, B, E, umin, umax, uPrevMin, uPrevMax)
    S.minHRep();
    H = S.A;
    g = S.b;
    margin = max(H*E*uPrevMin, H*E*uPrevMax);

    jointA = [H*A, H*B; 0, 0, 1; 0, 0, -1];
    jointb = [g-margin; umax; -umin];
    jointSet = Polyhedron('A', jointA, 'b', jointb);
    Pre = jointSet.projection(1:2);
    Pre.minHRep();
end

function Pdiff = pontryaginDifferenceWithInputSet(S, E, uPrevMin, uPrevMax)
    S.minHRep();
    H = S.A;
    g = S.b;
    margin = max(H*E*uPrevMin, H*E*uPrevMax);
    Pdiff = Polyhedron('A', H, 'b', g-margin);
    Pdiff.minHRep();
    if Pdiff.isEmptySet()
        error('X_R minus D is empty.');
    end
end

function [x, uFollower] = simulateController( ...
    mode, x0, uPredecessor, A, B, E, Qx, Qu, P, Tset, robustnessSet, ...
    xmin, xmax, umin, umax, N)

    Tset.minHRep();
    robustnessSet.minHRep();
    Hterminal = Tset.A;
    hterminal = Tset.b;
    Hrobust = robustnessSet.A;
    hrobust = robustnessSet.b;

    Nsim = numel(uPredecessor);
    x = zeros(2,Nsim+1);
    uFollower = zeros(1,Nsim);
    x(:,1) = x0;
    ops = sdpsettings('solver','quadprog','verbose',0);

    for k = 1:Nsim
        xPred = sdpvar(2,N+1,'full');
        uPred = sdpvar(1,N,'full');
        constraints = [xPred(:,1) == x(:,k)];
        objective = 0;

        for j = 1:N
            constraints = [constraints, ...
                xPred(:,j+1) == A*xPred(:,j) + B*uPred(j), ...
                xmin <= xPred(:,j) <= xmax, ...
                umin <= uPred(j) <= umax];
            objective = objective + xPred(:,j)'*Qx*xPred(:,j) ...
                                  + uPred(j)'*Qu*uPred(j);
        end

        if strcmp(mode,'rmpc')
            constraints = [constraints, Hrobust*xPred(:,2) <= hrobust];
        elseif strcmp(mode,'dmpc')
            constraints = [constraints, ...
                Hrobust*(xPred(:,2) + E*uPredecessor(k)) <= hrobust];
        else
            error('Unknown control mode.');
        end

        constraints = [constraints, ...
            xmin <= xPred(:,N+1) <= xmax, ...
            Hterminal*xPred(:,N+1) <= hterminal];
        objective = objective + xPred(:,N+1)'*P*xPred(:,N+1);

        diagnostics = optimize(constraints, objective, ops);
        if diagnostics.problem ~= 0
            error('%s became infeasible at k=%d: %s', upper(mode), k, diagnostics.info);
        end

        uFollower(k) = value(uPred(1));
        x(:,k+1) = A*x(:,k) + B*uFollower(k) + E*uPredecessor(k);
    end
end

function text = logicalText(value)
    if value
        text = 'true';
    else
        text = 'false';
    end
end
