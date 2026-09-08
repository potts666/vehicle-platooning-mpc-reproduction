%% Step 4: disturbance set D^i and robust control invariant set X_R^i
% This implements the set iteration behind Definition 2.5 in the paper:
% S_{j+1} = X intersect {x | there exists u in U such that
%                         A*x + B*u + d belongs to S_j for every d in D}.

clear; clc; close all;
mpt_init;

%% Vehicle i = 2 follows predecessor vehicle i-1 = 1
Ts = 1;
h = 1;
ds = 4;

A = [1 Ts; 0 1];
B = [-h*Ts - Ts^2/2; -Ts];
E = [Ts^2/2; Ts];

xmin = [-ds; -15];
xmax = [120;  15];
umin = -5;                  % U^2: follower input bound
umax =  3;

% Paper's input scaling: U^1 = U^2*C, C = diag(0.9, 0.9).
C = diag([0.9 0.9]);
Upredecessor = [umin umax]*C;
uPrevMin = Upredecessor(1);
uPrevMax = Upredecessor(2);

X = Polyhedron('lb', xmin, 'ub', xmax);

% D = E*U^1 is a line segment in the 2D error-state plane.
D = Polyhedron('V', [(E*uPrevMin)'; (E*uPrevMax)']);

fprintf('Follower input set U^2 = [%.3g, %.3g]\n', umin, umax);
fprintf('Predecessor input set U^1 = [%.3g, %.3g]\n', ...
    uPrevMin, uPrevMax);
fprintf('Disturbance endpoints: [%.3g, %.3g]^T and [%.3g, %.3g]^T\n', ...
    E(1)*uPrevMin, E(2)*uPrevMin, E(1)*uPrevMax, E(2)*uPrevMax);

%% Iteration for the maximal robust control invariant set
Xrobust = X;
maxIterations = 100;

for j = 1:maxIterations
    predecessorSet = robustPredecessor( ...
        Xrobust, A, B, E, umin, umax, uPrevMin, uPrevMax);

    Xnext = X.intersect(predecessorSet);
    Xnext.minHRep();

    if Xnext.isEmptySet()
        error('The robust control invariant set became empty at iteration %d.', j);
    end

    if Xnext == Xrobust
        Xrobust = Xnext;
        fprintf('Robust control invariant set converged after %d iterations.\n', j);
        break
    end

    Xrobust = Xnext;

    if j == maxIterations
        error(['The set did not converge within %d iterations. ', ...
            'Do not use this finite approximation as an exact X_R yet.'], ...
            maxIterations);
    end
end
fprintf('Number of X_R^i inequalities: %d\n', size(Xrobust.A,1));

%% Plot X and its robust control invariant subset X_R
figure('Color', 'w', 'Name', 'Step 4 - robust control invariant set');
X.plot('color', [0.75 0.75 0.75], 'linewidth', 1.2);
hold on
Xrobust.plot('color', [0.49 0.18 0.56], 'linewidth', 2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
grid on
xlabel('e_p [m]')
ylabel('e_v [m/s]')
title('State constraint set X^i and robust control invariant set X_R^i')
legend('State constraint set X^i', 'Robust control invariant set X_R^i', ...
    'Origin', 'Location', 'best')

%% Save all ingredients for RMPC construction in the next step
save('robust_sets_step4.mat', 'A', 'B', 'E', 'X', 'Xrobust', 'D', ...
    'xmin', 'xmax', 'umin', 'umax', ...
    'uPrevMin', 'uPrevMax', 'C', 'Ts', 'h', 'ds');

%% Local function: robust one-step predecessor of a polytope
function Pre = robustPredecessor(S, A, B, E, umin, umax, uPrevMin, uPrevMax)
% If S = {z | H*z <= g}, require A*x+B*u+d in S for all d=E*u_prev.
% Since u_prev lies in an interval, the worst value of every row H_r*d
% occurs at one of the two input bounds.

    S.minHRep();
    H = S.A;
    g = S.b;

    HdAtLower = H*E*uPrevMin;
    HdAtUpper = H*E*uPrevMax;
    disturbanceMargin = max(HdAtLower, HdAtUpper);

    % Build a polytope in [x; u] and project it onto x.
    jointA = [H*A, H*B; ...
              0,   0,   1; ...
              0,   0,  -1];
    jointb = [g - disturbanceMargin; umax; -umin];

    stateInputSet = Polyhedron('A', jointA, 'b', jointb);
    Pre = stateInputSet.projection(1:2);
    Pre.minHRep();
end
