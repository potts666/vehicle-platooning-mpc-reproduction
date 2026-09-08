%% Step 2: compute the LQR terminal set T^i
% The terminal set contains states from which the terminal LQR controller
% can keep satisfying both the state and input constraints forever.

clear; clc; close all;
mpt_init;

%% Same model and constraints as Step 1
Ts = 1;
h = 1;
ds = 4;

A = [1 Ts; 0 1];
B = [-h*Ts - Ts^2/2; -Ts];

xmin = [-ds; -15];
xmax = [120;  15];
umin = -5;
umax =  3;

Qx = eye(2);
Qu = 1;

%% Terminal LQR controller: u = -K*x
[K, P] = dlqr(A, B, Qx, Qu);
Acl = A - B*K;

%% States that are immediately admissible under the LQR law
% They must satisfy x in X and umin <= -K*x <= umax.
Hadmissible = [ eye(2); ...
               -eye(2); ...
               -K; ...
                K];
hadmissible = [ xmax; ...
               -xmin; ...
                umax; ...
               -umin];

Xadmissible = Polyhedron('A', Hadmissible, 'b', hadmissible);
Xadmissible.minHRep();

%% Maximal positively invariant set for x+ = Acl*x
% T_{j+1} = T_j intersection {x | Acl*x belongs to T_j}.
Tset = Xadmissible;
maxIterations = 100;

for j = 1:maxIterations
    predecessor = Tset.invAffineMap(Acl);
    Tnext = Tset.intersect(predecessor);
    Tnext.minHRep();

    if Tnext == Tset
        fprintf('Terminal set converged after %d iterations.\n', j);
        Tset = Tnext;
        break
    end

    Tset = Tnext;

    if j == maxIterations
        error('Terminal set did not converge within %d iterations.', maxIterations);
    end
end

if Tset.isEmptySet()
    error('The computed terminal set is empty. Check the model and constraints.');
end

fprintf('LQR terminal feedback: u = -[%g  %g] x\n', K(1), K(2));
fprintf('Number of terminal-set inequalities: %d\n', size(Tset.A,1));

%% Verify the one-step invariance condition at terminal-set vertices
Tset.minVRep();
terminalVertices = Tset.V';
uAtVertices = -K*terminalVertices;
nextVertices = Acl*terminalVertices;

assert(all(uAtVertices >= umin-1e-8) && all(uAtVertices <= umax+1e-8), ...
    'An LQR input constraint is violated at a terminal-set vertex.');

for vertexIndex = 1:size(nextVertices,2)
    assert(Tset.contains(nextVertices(:,vertexIndex)), ...
        'A successor state leaves the terminal set.');
end

%% Visualize the terminal set inside the immediately admissible region
figure('Color', 'w', 'Name', 'Step 2 - LQR terminal set');
Xadmissible.plot('color', [0.75 0.75 0.75], 'linewidth', 1.2);
hold on
Tset.plot('color', [0 0.45 0.74], 'linewidth', 2);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
grid on
xlabel('e_p [m]')
ylabel('e_v [m/s]')
title('LQR admissible region and terminal set T^i')
legend('Immediately LQR-admissible region', 'Terminal set T^i', 'Origin', ...
    'Location', 'best')

%% Save the terminal ingredients for the next MPC step
save('terminal_set_step2.mat', 'A', 'B', 'K', 'P', 'Acl', ...
    'Tset', 'Xadmissible', 'xmin', 'xmax', 'umin', 'umax', ...
    'Qx', 'Qu', 'Ts', 'h', 'ds');
